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


def get_plant_importance(graph):
    g_plants = [node for node in graph.nodes if graph.nodes[node]['bipartite'] == 'plant']

    degree_count = {t[0]: t[1] for t in graph.degree}
    degree_centrality = nx.algorithms.bipartite.degree_centrality(graph, g_plants)
    betweenness = nx.algorithms.bipartite.betweenness_centrality(graph, g_plants)
    closeness = nx.algorithms.bipartite.closeness_centrality(graph, g_plants)
    pagerank = nx.algorithms.pagerank(graph, weight='weight')

    return degree_count, degree_centrality, betweenness, closeness, pagerank


if __name__ == "__main__":
    print("analysis moved to mainframe")
