import re
import os
import pathlib

import requests as req
import csv

from taxonomy_info import get_nodelist_kingdoms

"""
    Parses network using the mangal web API, and converts them to the IWDB format,
    with metadata posted at the start of the worksheet.
"""

MANGAL_URL = "https://mangal.io/api/v2/"
TAX_TEMPLATE = "[{0}]"
TAX_PATH = "/home/userors/academics/netlab/mangal/taxonomy"


def mangal_base_request(request_specifier):
    """
    The base function for requesting information from mangal.
    :param request_specifier: string specifier about the request.
    :return: mangal response as a json dictionary.
    """
    request = MANGAL_URL + request_specifier
    response = req.get(request)
    if not response.ok:
        raise ConnectionError("request '{0}' failed with error code {1}.".format(request, response.status_code))
    return response.json()


def mangal_request_by_id(data_type, mangal_id):
    """
    Requests a single element of some data type from mangal. Always returns a single element
    in one page.
    :param data_type: data type of the requested element (node, edge, network, etc.).
    :param mangal_id: id of the requested element.
    :return: mangal response as a json dictionary.
    """
    return mangal_base_request("{data_type}/{value}".format(data_type=data_type, value=mangal_id))


def mangal_request_by_query(data_type, query_parameter):
    """
    Requests all elements matching some query from mangal. Queries can return multiple
    results spanning across several pages, therefore this function should always be called
    from a generator function like query_iterator.
    :param data_type: data type of the query elements.
    :param query_parameter: a dictionary of query parameters.
    :return: mangal response as a json dictionary.
    """
    query_text = "{data_type}?{query}".format(data_type=data_type,
        query="&".join(["=".join([k, str(v)]) for k, v in query_parameter.items()])
    )
    return mangal_base_request(query_text)


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


def get_network_vertices(network_id):
    """
    gets all vertices of a network, finding their kingdom either from a file
    or from requests to taxonomic web services.
    :param network_id: network id of the requested nodes.
    :return: a list of nodes as dictionaries (not namedtuples for kingdom completion mutability)
    """
    raw_vertices = []
    for vertex_json in query_iterator("node", {"network_id": network_id}):
        raw_vertices.append(vertex_json)

    node_kingdoms = get_nodelist_kingdoms(raw_vertices)
    return {v['id']: {"id": v['id'], "name": v['original_name'],
            "kingdom": node_kingdoms[v['id']]} for v in raw_vertices}


def get_network_edges(network_id):
    """
    gets all edges (interactions in mangal) of a network. Checks for network integrity
    in the process. In case a network has more then one interaction between a pair of nodes,
    or if it has an interaction with no value (or value 0), the network is deemed defective.
    :param network_id: network id of the requested interactions.
    :return: a dictionary of edges by node pairs.
    """
    edge_dict = dict()
    for e in query_iterator("interaction", {"network_id": network_id}):
        e_from = e['node_from']
        e_to = e['node_to']
        existing_id = edge_dict.get((e_from, e_to), -1)

        assert existing_id != e['id'], \
            "Multiple edges from {0} to {1}, keys {2} and {3}".format(
                e_from, e_to, e['id'], existing_id)
        assert e['value'] != 0, "Edge from {} to {} has value 0".format(e_from, e_to)

        # TODO consider changing edge_dict to frozenset keys.
        edge_dict[(e['node_from'], e['node_to'])] = {'value': e['value'],
                                                     'edge_id': e['id']}
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


def complete_kingdoms_by_interactions(vertices, edges):
    """
    Leverages the plant-pollinator network's bipartite structure to complete the kingdom
    classification of nodes with no kingdom. In case the network is not bipartite, it is
    not a valid plant-pollinator network and an assertion error is raised.
    :param vertices: network vertices
    :param edges: network edges
    :return: the number of nodes still missing a kingdom after inferring by
    interactions.
    """
    network_adjacency = construct_adjacency_list(edges)
    unresolved_vertices = len(vertices)
    current_unresolved = sum(
        map(lambda val: 1 if not val['kingdom'] else 0, vertices.values())
    )

    while unresolved_vertices > current_unresolved > 0:
        unresolved_vertices = current_unresolved
        current_unresolved = 0

        for vertex in vertices.values():
            if vertex['kingdom']:
                continue
            vertex_neighbors = network_adjacency.get(vertex['id'], [])
            neighbor_kingdoms = {vertices[nid]['kingdom'] for nid in vertex_neighbors} - {''}
            assert len(neighbor_kingdoms) < 2, "Network is not bipartite for vertex {0}".format(
                vertex['id'])

            if not neighbor_kingdoms:
                current_unresolved += 1
            else:
                vertex['kingdom'] = 'Animalia' if 'Plantae' in neighbor_kingdoms else 'Plantae'

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


if __name__ == "__main__":
    for network in query_iterator('network', {}):
        if not is_pollination_network(network['description']):
            continue

        nid = network['id']
        print("working on network {0}: {1}".format(nid, network['description']))
        try:
            nvertices = get_network_vertices(nid)
            nedges = get_network_edges(nid)
            print("{0} missing vertices for network {1}".format(
                complete_kingdoms_by_interactions(nvertices, nedges), nid
            ))
            for v in nvertices.values():
                if not v['kingdom']:
                    print(v)
        except AssertionError as e:
            print('network {0} error: {1}'.format(network['id'], e))
        print("=================================================")
