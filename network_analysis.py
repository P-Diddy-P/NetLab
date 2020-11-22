from collections import namedtuple

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def dataframe_to_network(dataframe):
    net_graph = nx.Graph()

    net_graph.add_nodes_from(list(dataframe.columns), bipartite='pollinator')
    net_graph.add_nodes_from(list(dataframe.index), bipartite='plant')

    for col_name in dataframe:
        for row_name in dataframe.index:
            if dataframe[col_name][row_name] > 0:
                net_graph.add_edge(col_name, row_name,
                                   weight=dataframe[col_name][row_name])
    return net_graph


def plant_measures(graph):
    """
    Calculates centrality measures for all plants in bipartite graph
    :param graph: Bipartite networkx graph
    :return: list of dictionaries with names and centrality measures of each species
    """
    g_plants = [node for node in graph.nodes if graph.nodes[node]['bipartite'] == 'plant']
    analysis_results = dict()

    degree_count = {t[0]: t[1] for t in graph.degree}
    degree_centrality = nx.algorithms.bipartite.degree_centrality(graph, g_plants)
    betweenness = nx.algorithms.bipartite.betweenness_centrality(graph, g_plants)
    closeness = nx.algorithms.bipartite.closeness_centrality(graph, g_plants)
    pagerank = nx.algorithms.pagerank(graph, weight='weight')

    for species in degree_count.keys():
        analysis_results[species] = {
            'dn': degree_count[species],
            'dc': degree_centrality[species],
            'cc': closeness[species],
            'bc': betweenness[species],
            'pr': pagerank[species]
        }
    return analysis_results


def rank_plants(network_table, interest):
    network_graph = dataframe_to_network(network_table)
    measures = plant_measures(network_graph)
    interest_rank = {sp: dict() for sp in interest}

    for sort_key in ['dn', 'dc', 'cc', 'bc', 'pr']:
        me_sort = [(sp, measures[sp][sort_key]) for sp in measures.keys()]
        print(sorted(me_sort, key=lambda e: e[1]))
        # create a list of species names and sort them by the importance
        # measure. After that, it is possible to fix the normalized rank of
        # species of interest by adding their distance from previous species
        # and subtracting the mean distance (not sure if necessary, add as param).
        # consider first counting polys in all networks.


if __name__ == "__main__":
    print("analysis moved to mainframe")
