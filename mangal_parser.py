import re

import requests as req
import openpyxl

from taxonomy_info import get_nodelist_kingdoms

"""
    Parses network using the mangal web API, and converts them to the IWDB format,
    with metadata posted at the start of the worksheet.
"""

MANGAL_URL = "https://mangal.io/api/v2/"
NAME_TEMPLATE = "[{0}]-{1}"
SAVE_PATH = "/home/userors/academics/netlab/mangal"


def mangal_base_request(request_specifier):
    request = MANGAL_URL + request_specifier
    response = req.get(request)
    if not response.ok:
        raise ConnectionError("request '{0}' failed with error code {1}.".format(request, response.status_code))
    return response.json()


def mangal_request_by_id(data_type, mangal_id):
    return mangal_base_request("{data_type}/{value}".format(data_type=data_type, value=mangal_id))


def mangal_request_by_query(data_type, query_parameter):
    query_text = "{data_type}?{query}".format(data_type=data_type,
                    query="&".join(["=".join([k, str(v)]) for k, v in query_parameter.items()])
                )
    return mangal_base_request(query_text)


def query_iterator(data_type, query):
    iter_query = {k: v for k, v in query.items()}
    iter_query['page'] = 0
    request_page = mangal_request_by_query(data_type, iter_query)
    while request_page:
        for entry in request_page:
            yield entry
        iter_query['page'] = iter_query['page'] + 1
        request_page = mangal_request_by_query(data_type, iter_query)


def get_network_metadata(network_id):
    return mangal_request_by_id("network", network_id)


def get_network_vertices(network_id):
    raw_vertices = []
    for vertex_json in query_iterator("node", {"network_id": network_id}):
        raw_vertices.append(vertex_json)

    node_kingdoms = get_nodelist_kingdoms(raw_vertices)
    return {v['id']: {"id": v['id'], "name": v['original_name'],
            "kingdom": node_kingdoms[v['id']]} for v in raw_vertices}


def get_network_edges(network_id):
    edge_dict = dict()
    for e in query_iterator("interaction", {"network_id": network_id}):
        e_from = e['node_from']
        e_to = e['node_to']
        existing_id = edge_dict.get((e_from, e_to), -1)

        assert existing_id != e['id'], \
            "Multiple edges from {0} to {1}, keys {2} and {3}".format(
                e_from, e_to, e['id'], existing_id)
        assert e['value'] != 0, "Edge from {} to {} has value 0".format(e_from, e_to)

        edge_dict[(e['node_from'], e['node_to'])] = {'value': e['value'],
                                                     'edge_id': e['id']}
    return edge_dict


def construct_adjacency_list(adjacency_matrix):
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
    including_criteria = ['pol[li]{1,3}nat', 'flower.+visit', 'flower.+interaction']
    excluding_criteria = ['food.web', 'pelagic', 'predat']

    included = any([re.search(pat, network_description.lower()) is not None for pat in including_criteria])
    excluded = any([re.search(pat, network_description.lower()) is not None for pat in excluding_criteria])
    return included and not excluded


def network_to_csv(vertices, edges, network_id):
    network_metadata = get_network_metadata(network_id)
    network_filename = NAME_TEMPLATE.format(network_id, network_metadata['name'])

    # TODO rows=plants and columns=pollinators
    pollinator_headers = []
    plant_headers = []
    with open(SAVE_PATH + network_filename, 'w+', newline='') as net_csv:
        net_csv.writerow([len(edges)].extend(pollinator_headers))

        raise NotImplementedError


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
