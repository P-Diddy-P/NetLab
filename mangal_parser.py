import re

import requests as req
import openpyxl

from taxonomy_info import get_nodelist_kingdoms

"""
    Parses network using the mangal web API, and converts them to the IWDB format,
    with metadata posted at the start of the worksheet.
"""

MANGAL_URL = "https://mangal.io/api/v2/"
SAVE_PATH = "C:\\Personal\\University\\Lab\\Mangal\\{filename}.xlsx"


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
        iter_query['page'] = iter_query.get('page', 0) + 1
        request_page = mangal_request_by_query(data_type, iter_query)


def get_network_metadata(network_id):
    return mangal_request_by_id("network", network_id)


def get_network_vertices(network_id):
    raw_vertices = []
    for vertex_json in query_iterator("node", {"network_id": network_id}):
        raw_vertices.append(vertex_json)

    node_kingdoms = get_nodelist_kingdoms(raw_vertices)
    return [{"id": v['id'], "name": v['original_name'],
            "kingdom": node_kingdoms[v['id']]} for v in raw_vertices]


def get_network_edges(network_id):
    edge_dict = dict()
    for e in query_iterator("interaction", {"network_id": network_id}):
        assert (e['node_from'], e['node_to']) not in edge_dict, \
            "Multiple edges from {} to {}".format(e['node_from'], e['node_to'])
        assert e['value'] != 0, \
            "Edge from {} to {} has value 0".format(e['node_from'], e['node_to'])

        edge_dict[(e['node_from'], e['node_to'])] = e['value']
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
    current_unresolved = len({v['id'] if v['kingdom'] else '' for v in vertices} - {''})

    while unresolved_vertices > current_unresolved > 0:
        unresolved_vertices = current_unresolved
        current_unresolved = 0
        for vertex in vertices:
            if not vertex['kingdom']:
                continue
            vertex_neighbors = network_adjacency.get(vertex['id'], [])
            if not vertex_neighbors:  # in case of isolated nodes, see node 4497 in network 27
                continue
            neighbor_kingdoms = {neighbor['kingdom'] for neighbor in vertex_neighbors} - {""}
            assert len(neighbor_kingdoms) < 2, "Network is not bipartite!"
            if not neighbor_kingdoms:
                current_unresolved += 1
            else:
                vertex['kingdom'] = 'Animalia' if 'Plantae' in neighbor_kingdoms else 'Plantae'
                print("Vertex {0} given kingdom {1} by connections".format(vertex['original_name'],
                                                                           vertex['kingdom']))


def fill_table_metadata(sheet, metadata):
    sheet['A1'] = str(metadata)
    sheet['C1'] = "plant_ge"
    sheet['C2'] = "plant_sp"
    sheet['C3'] = "no."
    sheet['A3'] = "pol_ge"
    sheet['B3'] = "pol_sp"


def fill_table_vertices(sheet, vertex_list):
    sheet_index = 4  # note that ws.cell function is 1 indexed
    for v in vertex_list:
        sheet.cell(row=sheet_index, column=1, value=v['genus'])
        sheet.cell(row=sheet_index, column=2, value=v['species'])
        sheet.cell(row=sheet_index, column=3, value=v['id'])

        sheet.cell(row=1, column=sheet_index, value=v['genus'])
        sheet.cell(row=2, column=sheet_index, value=v['species'])
        sheet.cell(row=3, column=sheet_index, value=v['id'])

        sheet_index += 1


def fill_table_edges(sheet, edge_dict, row_length):
    for i in range(4, 4 + row_length):
        for j in range(4, 4 + row_length):
            pair_id = sheet.cell(row=i, column=3).value, sheet.cell(row=3, column=j).value
            sheet.cell(row=i, column=j, value=edge_dict.get(pair_id, 0))


def construct_network_excel(network_id, save_path=SAVE_PATH):
    try:
        network_metadata = get_network_metadata(network_id)
        network_vertices = get_network_vertices(network_id)
        network_edges = get_network_edges(network_id)
    except AssertionError:
        print("Error parsing network ", network_id)
        raise

    network_book = openpyxl.Workbook()
    network_sheet = network_book.active

    fill_table_metadata(network_sheet, network_metadata)
    fill_table_vertices(network_sheet, network_vertices)
    fill_table_edges(network_sheet, network_edges, len(network_vertices))
    network_book.save(save_path.format(filename=network_metadata['name']))


def is_pollination_network(network_description):
    including_criteria = ['pol[li]{1,3}nat', 'flower.+visit', 'flower.+interaction']
    excluding_criteria = ['food.web', 'pelagic', 'predat']

    included = any([re.search(pat, network_description.lower()) is not None for pat in including_criteria])
    excluded = any([re.search(pat, network_description.lower()) is not None for pat in excluding_criteria])
    return included and not excluded


if __name__ == "__main__":
    for v in get_network_vertices(27):
        print(v)
    print("======================================================")
    print("======================================================")
    for e, v in get_network_edges(27).items():
        print("{f} - {t} :: {value}".format(f=e[0], t=e[1], value=v))
