from sys import argv
import pathlib
import os

import pandas as pd
import numpy as np
import scipy.stats as stats
from matplotlib import pyplot as plt
from Bio import Phylo

from polyploid_species import PolyploidDictionary
from duplicate_finder import compare_all_networks
from network_analysis import rank_graph, network_nodf, permutation_test, weighted_permutation_score
from phylogenetic_tree import generate_sparse_tree, closest_leaf_distances


def unidentified_rate(net_table):
    """
    Counts the number of COMPLETELY unidentified species of plants and pollinators
    in a network and divides by the total species number. Assumes such species are
    named `unidentified` in the network, and are wholly unknown (i.e. the genus is
    also unknown).
    :param net_table: pandas table of a pollinator network
    :return: tuple of unknown plant and pollinator rates.
    """
    unidentified_plants = sum(1 if p.lower().startswith('unidentified') else 0
                              for p in net_table.index)
    unidentified_pols = sum(1 if p.lower().startswith('unidentified') else 0
                            for p in net_table.columns)
    return unidentified_plants / len(net_table.index), unidentified_pols / len(net_table.columns)


def clean_networks(networks, plant_threshold=0.25, pol_threshold=1.0,
                   clean_plants=False, clean_pollinators=False):
    partial_networks = set()

    for net_name in networks.keys():
        ref_net = networks[net_name]

        try:
            missing_plants, missing_pols = unidentified_rate(ref_net)
        except Exception:
            print(net_name, "\n", ref_net)
            raise

        if clean_pollinators or clean_plants:
            known_cols = [not (c.lower().startswith('unidentified') and clean_pollinators)
                          for c in ref_net.columns]
            known_rows = [not (r.lower().startswith('unidentified') and clean_plants)
                          for r in ref_net.index]

            dropped_net = ref_net.loc[known_rows, known_cols]
            networks[net_name] = dropped_net

        if missing_plants > plant_threshold or missing_pols > pol_threshold:
            partial_networks.add(net_name)

    return partial_networks


def iter_sources(source_directories):
    for src in source_directories:
        for file in os.listdir(src):
            if file.endswith('csv'):
                yield src.joinpath(file)


def extract_networks(source_directories):
    all_networks = dict()
    duplicate_names = []
    for net_file in iter_sources(source_directories):
        if net_file.stem not in all_networks:
            all_networks[net_file.stem] = pd.read_csv(net_file, index_col=0)
        else:
            duplicate_names.append(net_file)

    return all_networks, duplicate_names


def permute_table(table, row_order, column_order):
    ordered = table.reindex(row_order, axis='index')
    ordered = ordered.reindex(column_order, axis='columns')
    return ordered


def generate_network_sparse_trees(networks, phylo_tree, tree_path):
    for name, network in networks.items():
        try:
            plant_species = list(network.index)
            sparse_tree_path = tree_path.parent.joinpath(f"sparse/{name}")
            generate_sparse_tree(phylo_tree, plant_species, sparse_tree_path)
        except Exception as e:
            print(f"failed to create sparse tree for network {name} [ERROR: {e}]")
        else:
            print(f"creates sparse tree for network {name}")


def analyze_networks(networks, polyploids, net_tree_path, measure='pg', missing_tolerance=0.25):
    net_rank_base, net_rank_permutation = {}, {}
    net_nodf_base, net_nodf_permutation = {}, {}
    polyploid_distances, diploid_distances = [], []

    for net_name, net_table in networks.items():
        net_poly = polyploids[net_name]
        r_plants, r_pols = rank_graph(net_table, measure_by=measure)

        ordered_table = permute_table(
            net_table,
            sorted(r_plants.keys(), key=r_plants.get, reverse=True),
            sorted(r_pols, key=r_pols.get, reverse=True)
        )
        total_nodf, plant_nodf, pol_nodf = network_nodf(ordered_table, correction='sqrt')

        try:  # rank permutation tests
            rank_base, rank_permutations = permutation_test(r_plants, net_poly)
            net_rank_base[net_name] = rank_base
            net_rank_permutation[net_name] = rank_permutations
        except FloatingPointError:
            print(f"error calculating rank permutations in network {net_name}, continuing...")

        try:  # nodf permutation tests
            nodf_base, nodf_permutation = permutation_test(plant_nodf, net_poly)
            net_nodf_base[net_name] = nodf_base
            net_nodf_permutation[net_name] = nodf_permutation
        except FloatingPointError:
            print(f"error calculating nodf permutations in network {net_name}, continuing...")

        try:  # phylogenetic tree tests
            sparse_tree = Phylo.read(net_tree_path.joinpath(net_name, "tree.tre"), 'newick')
            with open(net_tree_path.joinpath(net_name, "missing"), mode='r') as fd:
                missing_species = fd.read().split(',')
                missing_rate = len(missing_species) / len(net_table.index)
                if missing_rate > missing_tolerance:
                    print(f"network {net_name} has too many missing species for accurate distance measures.")
                    continue

            # Phylo.draw(sparse_tree)
            for clade, distance in closest_leaf_distances(sparse_tree).items():
                species = clade.name.replace('_', ' ')
                if species in net_poly:
                    polyploid_distances.append(distance)
                else:
                    diploid_distances.append(distance)
        except (FloatingPointError, FileNotFoundError):
            print(f"Error getting closest leaf distances for network: {net_name}")

    return net_rank_base, net_rank_permutation, net_nodf_base, net_nodf_permutation, \
        polyploid_distances, diploid_distances


def visualize_results(networks, polyploids, rank_base, rank_perm, nodf_base, nodf_perm,
                      distance_poly, distance_di, result_path):
    rank_over_base = {
        net_name: (sum([1.0 for perm in permutations if perm > rank_base[net_name]]) / 10000.0) for
        net_name, permutations in rank_perm.items()
    }
    stats.probplot(list(rank_over_base.values()), dist='uniform', plot=plt)
    plt.show()

    rank_mean_overbase_no_fix = sum([
        1.0 for i in range(10000) if weighted_permutation_score(rank_base, rank_perm, i, networks)
    ]) / 10000.0
    rank_mean_overbase_log_fix = sum([
        1.0 for i in range(10000) if weighted_permutation_score(rank_base, rank_perm, i, networks, weight=np.log)
    ]) / 10000.0
    rank_mean_overbase_line_fix = sum([
        1.0 for i in range(10000) if weighted_permutation_score(rank_base, rank_perm, i, networks, weight='id')
    ]) / 10000.0
    print(f"rank mean overbase no fix {rank_mean_overbase_no_fix}")
    print(f"rank mean overbase log fix {rank_mean_overbase_log_fix}")
    print(f"rank mean overbase linear fix {rank_mean_overbase_line_fix}")

    nodf_over_base = {
        net_name: (sum([1.0 for perm in permutations if perm > nodf_base[net_name]]) / 10000.0) for
        net_name, permutations in nodf_perm.items()
    }
    stats.probplot(list(nodf_over_base.values()), dist='uniform', plot=plt)
    plt.show()

    nodf_mean_overbase_no_fix = sum([
        1.0 for i in range(10000) if weighted_permutation_score(nodf_base, nodf_perm, i, networks)
    ]) / 10000.0
    nodf_mean_overbase_log_fix = sum([
        1.0 for i in range(10000) if weighted_permutation_score(nodf_base, nodf_perm, i, networks, weight=np.log)
    ]) / 10000.0
    nodf_mean_overbase_line_fix = sum([
        1.0 for i in range(10000) if weighted_permutation_score(nodf_base, nodf_perm, i, networks, weight='id')
    ]) / 10000.0
    print(f"nodf mean overbase no fix {nodf_mean_overbase_no_fix}")
    print(f"nodf mean overbase log fix {nodf_mean_overbase_log_fix}")
    print(f"nodf mean overbase linear fix {nodf_mean_overbase_line_fix}")

    # size_rank = [(len(net_table.index), rank_over_base[net_name]) for net_name, net_table in networks.items()]
    # plt.scatter([e[0] for e in size_rank], [e[1] for e in size_rank])
    # plt.show()

    # after that, perform two sample KS test on polyploid vs diploids shortest distances (consider
    # both poly<>di and poly>=di hypotheses).
    print(len(distance_poly), sorted(distance_poly))
    print(len(distance_di), sorted(distance_di))
    print(stats.kstest(distance_poly, distance_di, alternative='two-sided'))
    print(stats.kstest(distance_poly, distance_di, alternative='greater'))
    print(stats.kstest(distance_poly, distance_di, alternative='less'))

    # finally, perform linear regression of shortest distances to max divided importance (requires
    # calculating importance again or moving from analyze_networks


if __name__ == "__main__":
    np.seterr(all='raise')

    analysis_results_path = pathlib.Path(argv[1])
    polyploid_path = pathlib.Path(argv[2])
    tree_path = pathlib.Path(argv[3]).joinpath('ALLMB.tre')
    network_paths = [pathlib.Path(path_str) for path_str in argv[4:]]

    polydict = PolyploidDictionary(polyploid_path)
    phylo_tree = Phylo.read(tree_path, format='newick')

    network_tables, duplicate_network_names = extract_networks(network_paths)
    invalid_networks = clean_networks(network_tables)

    duplicates = compare_all_networks(network_tables, 0.05, drop_early=invalid_networks)
    print(f"{len(duplicates)} duplicates found")
    network_polyploids = dict()
    for name, table in network_tables.items():
        polies = {pn for pn in table.index if polydict.test_ploidy(pn)}
        if len(polies) > 0:
            network_polyploids[name] = polies

    networks_to_analyze = {name: table for name, table in network_tables.items() if
                           name in network_polyploids and
                           name not in duplicates and
                           name not in invalid_networks}

    """ SPARSE TREE CREATION IS TIME CONSUMING AND ALREADY DONE
    generate_network_sparse_trees(networks_to_analyze, phylo_tree, tree_path)
    """

    print(f"analyzing {len(networks_to_analyze)} networks")
    rank_b, rank_p, nodf_b, nodf_p, dist_p, dist_d = analyze_networks(
        networks_to_analyze, network_polyploids, tree_path.parent.joinpath("sparse"), missing_tolerance=0.25
    )

    visualize_results(networks_to_analyze, network_polyploids, rank_b, rank_p, nodf_b, nodf_p,
                      dist_p, dist_d, analysis_results_path)
