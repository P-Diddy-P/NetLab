from sys import argv
import re
import os
import pathlib
import csv

import requests as req

from taxonomy_info import get_nodelist_kingdoms

"""
    Parses network using the mangal web API, and converts them to the IWDB format,
    with metadata posted at the start of the worksheet.
"""

MANGAL_URL = "https://mangal.io/api/v2/"


def mangal_base_request(request_specifier, session):
    """
    The base function for requesting information from mangal.
    :param request_specifier: string specifier about the request.
    :param session: a requests.session object to send requests from.
    :return: mangal response as a json dictionary.
    """
    request = MANGAL_URL + request_specifier
    if session:
        response = session.get(request)
    else:
        response = req.get(request)

    if not response.ok:
        raise ConnectionError("request '{0}' failed with error code {1}.".format(
            request, response.status_code))

    return response.json()


def mangal_request_by_id(data_type, mangal_id, session=None):
    """
    Requests a single element of some data type from mangal. Always returns a single element
    in one page.
    :param data_type: data type of the requested element (node, edge, network, etc.).
    :param mangal_id: id of the requested element.
    :param session: a requests.session object to send requests from.
    :return: mangal response as a json dictionary.
    """
    return mangal_base_request(
        "{data_type}/{value}".format(data_type=data_type, value=mangal_id), session
    )


def mangal_request_by_query(data_type, query_parameter, session=None):
    """
    Requests all elements matching some query from mangal. Queries can return multiple
    results spanning across several pages, therefore this function should always be called
    from a generator function like query_iterator.
    :param data_type: data type of the query elements.
    :param query_parameter: a dictionary of query parameters.
    :param session: a requests.session object to send requests from.
    :return: mangal response as a json dictionary.
    """
    query_text = "{data_type}?{query}".format(data_type=data_type,
        query="&".join(["=".join([key, str(val)]) for key, val in query_parameter.items()])
    )
    return mangal_base_request(query_text, session)


def query_iterator(data_type, query):
    """
    A generator function for mangal queries. Iterates over all entries from
    the query by the same ordering as they're retrieved from mangal, no matter
    over how many pages the query spans.
    :param data_type: data type of the query elements.
    :param query: a dictionary of query key-value elements (is not kwargs for a reason).
    :return: an iterator of query results.
    """
    iter_query = dict(query)
    iter_query['page'] = 0
    request_page = mangal_request_by_query(data_type, iter_query)
    while request_page:
        for entry in request_page:
            yield entry
        iter_query['page'] = iter_query['page'] + 1
        request_page = mangal_request_by_query(data_type, iter_query)


def get_network_metadata(network_id):
    """
    :param network_id: id of the requested network
    :return: network metadata such as description, geo-location, authors, etc.
    """
    return mangal_request_by_id("network", network_id)


def get_network_nodes(network_id, node_kingdoms=None, force_web=False):
    """
    gets all vertices of a network, finding their kingdom either from a file
    or from requests to taxonomic web services.
    :param network_id: network id of the requested nodes.
    :param node_kingdoms: a preloaded dictionary of node kingdoms. If no such dictionary exists,
    a new one will be generated from taxonomy web services.
    :param force_web: set to true to force using web services for node kingdoms.
    :return: a list of nodes as dictionaries (not namedtuples for kingdom completion mutability)
    """
    raw_nodes = []
    node_counter = 0
    for node_json in query_iterator("node", {"network_id": network_id}):
        raw_nodes.append(node_json)
        node_counter += 1

    if force_web or not node_kingdoms:
        node_kingdoms = get_nodelist_kingdoms(raw_nodes)
    all_nodes = {node['id']: {"id": node['id'], "name": node['original_name'],
                              "kingdom": node_kingdoms.get(node['id'], '')}
                 for node in raw_nodes}
    assert len(all_nodes) == node_counter, \
        f"total nodes should be {node_counter}, are {len(all_nodes)}"
    return all_nodes


def get_network_edges(network_id):
    """
    gets all edges (interactions in mangal) of a network. Checks for network integrity
    in the process. In case a network has more then one interaction between a pair of nodes,
    or if it has an interaction with no value (or value 0), the network is deemed defective.
    :param network_id: network id of the requested interactions.
    :return: a dictionary of edges by node pairs.
    """
    edge_dict = dict()
    edge_counter = 0
    for edge in query_iterator("interaction", {"network_id": network_id}):
        edge_counter += 1
        e_from = edge['node_to']
        e_to = edge['node_from']
        existing_id = edge_dict.get(frozenset((e_from, e_to)), -1)

        assert existing_id != edge['id'], \
            "Multiple edges from {0} to {1}, keys {2} and {3}".format(
                e_from, e_to, edge['id'], existing_id)
        assert edge['value'] != 0, "Edge from {} to {} has value 0".format(e_from, e_to)

        edge_dict[frozenset((edge['node_from'], edge['node_to']))] = \
            {'value': edge['value'], 'edge_id': edge['id']}  # , 'lid': edge_counter}

    assert edge_counter == len(edge_dict), f"total interactions should be {edge_counter}, are {len(edge_dict)}"
    return edge_dict


def construct_adjacency_list(adjacency_matrix):
    """
    construct an adjacency list representation of network interactions
    from an adjacency matrix, for easier searching of neighbours.
    :param adjacency_matrix: an adjacency matrix representation of a network.
    :return: a dictionary mapping between each node id and its neighbours.
    """
    network_adjacency = dict()
    for src, dest in adjacency_matrix.keys():
        if src in network_adjacency:
            network_adjacency[src].append(dest)
        else:
            network_adjacency[src] = [dest]

        if dest in network_adjacency:
            network_adjacency[dest].append(src)
        else:
            network_adjacency[dest] = [src]

    return network_adjacency


def complete_kingdoms_by_interactions(nodes, edges):
    """
    Leverages the plant-pollinator network's bipartite structure to complete the kingdom
    classification of nodes with no kingdom. In case the network is not bipartite, it is
    not a valid plant-pollinator network and an assertion error is raised.
    :param nodes: network vertices
    :param edges: network edges
    :return: the number of nodes still missing a kingdom after inferring by
    interactions.
    """
    network_adjacency = construct_adjacency_list(edges)
    unresolved_nodes = len(nodes)
    current_unresolved = sum(
        map(lambda val: 1 if not val['kingdom'] else 0, nodes.values())
    )

    while unresolved_nodes > current_unresolved > 0:
        unresolved_nodes = current_unresolved
        current_unresolved = 0

        for node in nodes.values():
            if node['kingdom']:
                continue
            node_neighbors = network_adjacency.get(node['id'], [])
            neighbor_kingdoms = {nodes[neighbor]['kingdom'] for neighbor in node_neighbors} - {''}
            assert len(neighbor_kingdoms) < 2, \
                "Network is not bipartite for vertex {0}: {1}".format(
                node['id'], {n: nodes[n]['kingdom'] for n in node_neighbors})

            if not neighbor_kingdoms:
                current_unresolved += 1
            else:
                node['kingdom'] = 'Animalia' if 'Plantae' in neighbor_kingdoms else 'Plantae'

    return current_unresolved


def is_pollination_network(network_description):
    """
    Superficial check for plant-pollinator networks according to the network's mangal
    description. The test is based on some patterns perceived from manually looking at
    network names, and could possibly be extended.
    :param network_description: description of the network in question.
    :return: true if network is a plant pollinator network, false otherwise.
    """
    including_criteria = ['pol[li]{1,3}nat', 'flower.+visit', 'flower.+interaction']
    excluding_criteria = ['food.web', 'pelagic', 'predat']

    included = any([re.search(pat, network_description.lower()) is not None for pat in including_criteria])
    excluded = any([re.search(pat, network_description.lower()) is not None for pat in excluding_criteria])
    return included and not excluded


def load_network_taxonomy(network_id, load_path):
    """
    Loads a network's taxonomy data from file (if one exists).
    :param network_id: network id of the loaded nodes.
    :param load_path: directory path to load from.
    :return: a dictionary of node ids mapped to kingdom if a
    taxonomy file exists, otherwise None.
    """
    taxonomy_filename = 'net{0}'.format(network_id)
    if taxonomy_filename not in os.listdir(load_path):
        return

    taxonomy_dict = dict()
    with open(load_path.joinpath(taxonomy_filename), 'r') as tax_file:
        plant_line = tax_file.readline()[:-1]
        for node_id in plant_line.split('|'):
            taxonomy_dict[int(node_id)] = 'Plantae'

        pollinator_line = tax_file.readline()[:-1]
        for node_id in pollinator_line.split('|'):
            taxonomy_dict[int(node_id)] = 'Animalia'

        unknown_line = tax_file.readline()[:-1]
        for node_id in unknown_line.split('|'):
            if node_id:
                taxonomy_dict[int(node_id)] = ''
    return taxonomy_dict


def store_network_taxonomy(network_id, network_nodes, store_path):
    """
    stores a network's taxonomy data in a file as node ids separated
    by pipe (|) characters. First row is for plants, second for pollinators
    and third for nodes with unknown taxonomy. After storing, the taxonomy
    data can be completed manually.
    :param network_id: network id of the stored vertices.
    :param store_path: directory path to store to.
    :param network_nodes: a list of vertices with taxonomy data
    """
    plants, pollinators, unknowns = [], [], []
    for node in network_nodes.values():
        store_id = str(node['id'])
        if node['kingdom'] == 'Plantae':
            plants.append(store_id)
        elif node['kingdom'] == 'Animalia':
            pollinators.append(store_id)
        else:  # no known kingdom for node.
            unknowns.append(store_id)

    store_path = pathlib.Path(store_path).joinpath('net{0}'.format(network_id))
    with open(store_path, 'w+') as store_file:
        store_file.write('|'.join(plants) + '\n')
        store_file.write('|'.join(pollinators) + '\n')
        store_file.write('|'.join(unknowns) + '\n')


def network_to_csv(nodes, edges, csv_path):
    """
    converts the network into a pandas friendly csv format. The format conforms
    to WoL csv format with plants as rows starts and pollinators as column heads
    :param nodes: node list of the selected network.
    :param edges: interaction dict of the selected network.
    :param csv_path: path to save the csv.
    """
    plant_nodes = [nodes[nid] for nid in nodes if nodes[nid]['kingdom'] == 'Plantae']
    pollinator_nodes = [nodes[nid] for nid in nodes if nodes[nid]['kingdom'] == 'Animalia']
    with open(csv_path, 'w+', newline='') as csvfile:
        net_writer = csv.writer(csvfile, delimiter=',')
        net_writer.writerow([''] + [poll['name'] for poll in pollinator_nodes])
        for plant in plant_nodes:
            plant_row = [plant['name']]
            for poll in pollinator_nodes:
                node_pair = frozenset((plant['id'], poll['id']))
                edge_value = 0
                if node_pair in edges:
                    edge_value = edges[node_pair]['value']
                plant_row.append(edge_value)
            net_writer.writerow(plant_row)


def construct_network(network_id, base_dir, force_web=False, force_csv=False):
    """
    Perform the whole network construction pipeline, getting nodes, taxonomy
    and edges, then saving the network to csv.
    :param network_id: mangal id of the requested network.
    :param base_dir: base dir for storing network csv or taxonomy files from.
    :param force_web: get node taxonomy from web even if a taxonomy file exists.
    :param force_csv: create a csv even if one already exists (overwrite existing csv).
    :return: the network's nodes and edges.
    """
    net_path = pathlib.Path(base_dir).joinpath('netcsv', '{0}.csv'.format(
        get_network_metadata(network_id)['name']))
    if not force_csv and os.path.isfile(net_path):
        print("net csv found")  # TODO remove
        return net_path  # currently construct_network return values are not used, mostly symbolic

    preloaded_kingdoms = load_network_taxonomy(network_id, base_dir.joinpath('taxonomy'))
    net_nodes = get_network_nodes(
        network_id, node_kingdoms=preloaded_kingdoms, force_web=force_web
    )
    net_edges = get_network_edges(network_id)
    if force_web or not preloaded_kingdoms:
        missing_kingdoms = complete_kingdoms_by_interactions(net_nodes, net_edges)
        store_network_taxonomy(network_id, net_nodes, base_dir.joinpath('taxonomy'))
    else:
        missing_kingdoms = sum(
            map(lambda val: 1 if not val['kingdom'] else 0, net_nodes.values())
        )

    assert missing_kingdoms < 1, "network {0} has {1} nodes with no kingdom." \
                                 "Can't write to csv".format(network_id, missing_kingdoms)

    network_to_csv(net_nodes, net_edges, net_path)
    return net_nodes, net_edges


if __name__ == "__main__":
    project_path = pathlib.Path(argv[1])
    for i in range(1, 21):
        print("=" * 100)
        print(f"STARTING NETWORK RUN {i}")
        print("=" * 100)

        failed_nets = 0

        for network in query_iterator('network', {}):
            if not is_pollination_network(network['description']):
                continue

            nid = network['id']
            print("working on network {0}: {1}".format(nid, network['description']))
            try:
                construct_network(nid, project_path, force_web=True)
            except AssertionError as e:
                print(">>>>>>>>>>>>>>>>>>")
                print('network {0} error: {1}'.format(network['id'], e))
                print("<<<<<<<<<<<<<<<<<<")
                failed_nets += 1
            print("-" * 100)

        print(f"run {i} done with {failed_nets} failures.")
