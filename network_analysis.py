from itertools import chain, repeat

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def dataframe_to_network(dataframe, visualize=False):
    net_graph = nx.Graph()

    net_graph.add_nodes_from(list(dataframe.columns), bipartite='pollinator')
    net_graph.add_nodes_from(list(dataframe.index), bipartite='plant')

    for col_name in dataframe:
        for row_name in dataframe.index:
            if dataframe[col_name][row_name] > 0:
                strength = dataframe[col_name][row_name]
                net_graph.add_edge(col_name, row_name, weight=strength)

    if visualize:
        nx.draw(net_graph, with_labels=True)
        plt.show()
    return net_graph


def musrank_implementation(table, n_iterations=100, delta=10e-6):
    # fails to reach adequate convergence properties, dropped.
    """
    Implements a simplified version of the MusRank algorithm
    as seen in (DOI: 10.1038/srep08182).
    :param table: pandas dataframe of network.
    :param n_iterations: number of iterations to run MusRank.
    :param delta: threshold for accepting steady state.
    :return: a dictionary assigning importance to each node in the network.
    """
    no_isolated = table[(table.T != 0).any()]  # remove isolated nodes
    no_isolated = no_isolated.loc[:, (no_isolated != 0).any()]
    isolated_plants, isolated_pollinators = set(table.index) - set(no_isolated.index), \
                                            set(table.columns) - set(no_isolated.columns)

    n_plants, n_pols = no_isolated.shape
    zo_table = no_isolated.astype(bool).astype(float).to_numpy()
    plant_v, pol_i = np.ones(n_plants), np.ones(n_pols)

    for i in range(n_iterations):
        next_i = zo_table.transpose().dot(plant_v)
        normed_i = next_i / next_i.mean()

        inverse_i = np.divide(1, pol_i, out=np.ones_like(pol_i), where=(pol_i != 0))
        next_v = 1 / (zo_table.dot(inverse_i))
        normed_v = next_v / next_v.mean()

        iter_delta = np.linalg.norm(normed_v) + np.linalg.norm(normed_i)
        plant_v, pol_i = normed_v, normed_i
        if iter_delta < delta * n_pols * n_plants:
            print(f"Convergence reached in iteration {i}")
            break

    plant_rank = zip(no_isolated.index, plant_v)
    pol_rank = zip(no_isolated.columns, pol_i)
    plant_zeros = zip(isolated_plants, repeat(0))
    pol_zeros = zip(isolated_pollinators, repeat(0))

    net_rank = {k: v for k, v in chain.from_iterable([pol_rank, plant_rank, plant_zeros, pol_zeros])}
    return net_rank


def calculate_measure(graph, plant_side, measure='pg', weight='weight', table=None):
    if measure == 'dn':
        return {t[0]: t[1] for t in graph.degree}
    if measure == 'dc':
        return nx.algorithms.bipartite.degree_centrality(graph, plant_side)
    if measure == 'cc':
        return nx.algorithms.bipartite.closeness_centrality(graph, plant_side)
    if measure == 'bc':
        return nx.algorithms.bipartite.betweenness_centrality(graph, plant_side)
    if measure == 'pg':
        return nx.algorithms.pagerank(graph, weight=weight)
    if measure == 'mr' and table is not None:
        return musrank_implementation(table)

    raise ValueError('Unknown measure: {}'.format(measure))


def rank_graph(network_table, measure_by):
    """
    Calculates a centrality measure for all nodes in a network, and sorts them by that measure,
    returns a subset of the nodes if exist.
    :param network_table: pandas table to convert into networkx graph.
    :param measure_by: contrality measure, select one from calculate_measure if/else
    :return: two dictionaries of (species, rank), one for plants, one for pollinators
    """
    network_graph = dataframe_to_network(network_table)
    network_plants = set(network_table.index)
    network_pollinators = set(network_table.columns)
    ranked_graph = calculate_measure(network_graph, network_plants, measure_by,
                                     table=network_table)

    ranked_plants = dict((k, ranked_graph[k]) for k in network_plants)
    ranked_pollinators = dict((k, ranked_graph[k]) for k in network_pollinators)
    return get_fractional_indices(ranked_plants), get_fractional_indices(ranked_pollinators)


def get_fractional_indices(ranks, species_to_get=None):
    """
    Takes a ranks dictionary of species from a network, normalizes them and returns
    a subset.
    :Note this function is deprecated, doesn't effect log-odds measures.
    :param ranks: a dictionary of species and their rankings.
    :param species_to_get: which species to return, set to None to return all.
    :return: a dictionary of normalized measures for species_to_get.
    """
    if not species_to_get:
        species_to_get = ranks.keys()

    max_importance = max(ranks.values())
    return {k: v / max_importance for k, v in ranks.items() if k in species_to_get}


def network_nodf(network_table, correction=lambda n: n):
    """
    Calculate network NODF index for the entire network, as well as
    per species NODF contributions.

    Base algorithm: http://guimaraes.bio.br/032.pdf
    Reference implementation: https://github.com/tsakim/nestedness/blob/master/nestedness_calculator.py
    :param network_table: a pandas dataframe of the network
    :return: A 2-tuple of (network NODF score, species relative NODF dictionary)
    """
    if correction == 'sqrt':
        correction = lambda n: np.sqrt(n)
    elif correction is None:
        correction = lambda n: 1

    zo_table = network_table.astype(bool).astype(float).to_numpy()

    # Calculate row N_paired
    row_overlap = np.dot(zo_table, zo_table.T)
    row_degrees = zo_table.sum(axis=1)
    row_degrees_matrix = row_degrees * np.ones_like(row_overlap)

    df_paired = (row_degrees > row_degrees[:, np.newaxis])
    valid_matrix = np.multiply(df_paired, np.tril(np.ones_like(row_overlap), k=-1))
    row_partial_overlap = np.divide(row_overlap, row_degrees_matrix,
                                    out=np.zeros_like(row_overlap),
                                    where=(row_degrees_matrix != 0))

    row_nodf = np.dot(row_partial_overlap, valid_matrix).diagonal()

    # Calculate column N_paired
    col_overlap = np.dot(zo_table.T, zo_table)
    col_degrees = zo_table.sum(axis=0)
    col_degrees_matrix = col_degrees * np.ones_like(col_overlap)

    df_paired_col = (col_degrees > col_degrees[:, np.newaxis])
    valid_matrix = np.multiply(df_paired_col, np.tril(np.ones_like(col_overlap), k=-1))
    col_partial_overlap = np.divide(col_overlap, col_degrees_matrix,
                                    out=np.zeros_like(col_overlap),
                                    where=(col_degrees_matrix != 0))

    col_nodf = np.dot(col_partial_overlap, valid_matrix).diagonal()

    # Calculate final NODF and create dict with relative NODF importance
    table_shape = np.array(zo_table.shape)
    nodf_final = (col_nodf.sum() + row_nodf.sum()) / (np.dot(table_shape, table_shape - 1) / 2)
    plant_nodf = {network_table.index[i]: row_nodf[i] / correction(len(network_table.index) - i - 1)
                  for i in range(len(network_table.index) - 1)}
    pol_nodf = {network_table.columns[i]: col_nodf[i] / correction(len(network_table.columns) - i - 1)
                for i in range(len(network_table.columns) - 1)}

    plant_nodf[network_table.index[-1]] = 0
    pol_nodf[network_table.columns[-1]] = 0
    return nodf_final, plant_nodf, pol_nodf


def permutation_test(species_measures, base_polyploids, permutations=10000, tolerance=100):
    """
    Calculates the log odds between the mean measure of diploids and mean measure
    of polyploids, this is the base result. Then permutes the identity of polyploids
    ;permutations; times and calculates those simulated results.
    :param species_measures: a dictionary containing each plant species in
    the network and its measure.
    :param base_polyploids: a set of the network's polyploid species.
    :param permutations: number of permutations to perform. Defaults to 10000
    :param tolerance: number of floating point errors tolerated before raising error.
    :return: a 2-tuple containing the base result and a list of length ;permutation;
    with simulated results.
    """
    all_species, n_polyploids = list(species_measures.keys()), len(base_polyploids)

    base_diploids = set(all_species) - base_polyploids
    base_polyploid_mean = np.mean([species_measures[sp] for sp in base_polyploids])
    base_diploid_mean = np.mean([species_measures[sp] for sp in base_diploids])
    base_statistic = np.log(base_polyploid_mean / base_diploid_mean)

    permutation_statistics = []
    error_count = 0

    while len(permutation_statistics) < permutations:
        p_polyploids = set(np.random.choice(all_species, n_polyploids, replace=False))
        p_diploids = set(all_species) - p_polyploids
        try:
            p_polyploid_mean = np.mean([species_measures[sp] for sp in p_polyploids])
            p_diploid_mean = np.mean([species_measures[sp] for sp in p_diploids])

            permutation_statistics.append(np.log(p_polyploid_mean / p_diploid_mean))
        except FloatingPointError:
            error_count += 1
            if error_count >= tolerance:
                raise

    return base_statistic, permutation_statistics


def weighted_permutation_score(base, permutations, index, tables, weight=lambda n: 1):
    if weight == 'id':
        weight = lambda n: n

    weighted_base_mean = sum([
        net_base_score * weight(len(tables[net_name].index)) for net_name, net_base_score
        in base.items()
    ])
    weighted_perm_mean = sum([
        net_permutations[index] * weight(len(tables[net_name].index)) for net_name, net_permutations
        in permutations.items()
    ])

    return weighted_perm_mean > weighted_base_mean


if __name__ == "__main__":
    example = {
        'a': 100,
        'c': 80,
        'e': 55,
        'b': 122,
        'd': 43,
        'f': 500,
        'g': 32
    }
    print(example)
    print(get_fractional_indices(example, ['g']))
