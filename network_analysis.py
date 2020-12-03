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
                strength = dataframe[col_name][row_name]
                net_graph.add_edge(col_name, row_name, weight=strength)
    return net_graph


def calculate_measure(graph, plant_side, measure='pg', weight='weight'):
    print([node for node in graph.nodes if node in plant_side])
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

    raise ValueError('Unknown measure: {}'.format(measure))


def rank_plants(network_table, measure_by):
    """
    Calculates a centrality measure for all nodes in a network, and sorts them by that measure,
    returns a subset of the nodes if exist.
    :param network_table: pandas table to convert into networkx graph.
    :param measure_by: contrality measure, select one from calculate_measure if/else
    :return: a sorted list containing pairs of (plant species, centrality value)
    """
    network_graph = dataframe_to_network(network_table)
    network_plants = set(network_table.index)
    ranked_graph = calculate_measure(network_graph, network_plants, measure_by)

    ranked_plants = [(k, v) for k, v in ranked_graph.items() if k in network_plants]
    ranked_plants.sort(key=lambda elem: elem[1])
    return ranked_plants


if __name__ == "__main__":
    print("analysis moved to mainframe")
