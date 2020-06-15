import requests as req
import openpyxl

from taxonomy_info import get_kingdom

"""
    Parses network using the mangal web API, and converts them to the IWDB format,
    with metadata posted at the start of the worksheet.
"""

# TODO not all species in mangal are named and identified correctly, so in order to
# TODO coerce IWDB format, all unidentified species should be derived from their
# TODO connections

BY_ID_URL = "https://mangal.io/api/v2/{data_type}/{value}"
QUERY_URL = "https://mangal.io/api/v2/{data_type}?{query_by}={value}"
SAVE_PATH = "C:\\Personal\\University\\Lab\\Mangal\\{filename}.xlsx"


def mangal_request(data_type, query_parameter, request_value, is_query=True):
    if is_query:
        request = QUERY_URL.format(data_type=data_type, query_by=query_parameter, value=request_value)
    else:
        request = BY_ID_URL.format(data_type=data_type, value=request_value)
    response = req.get(request)

    if response.status_code != 200:
        raise ConnectionError("request '{0}' failed with error code {1}.".format(request, response.status_code))
    return response.json()


def get_network_metadata(network_id):
    # TODO consider slicing only part of the metadata (id, name, desc. etc.)
    return mangal_request("network", None, network_id, False)


def get_network_vertices(network_id):
    vertex_json = mangal_request("node", "network_id", network_id)
    return [{"id": v['id'], "genus": v['original_name'].split(' ')[0],
             "species": v['original_name'].split(' ')[1]} for v in vertex_json]


def get_network_edges(network_id):
    edge_json = mangal_request("interaction", "network_id", network_id)
    edge_dict = dict()
    for e in edge_json:
        if (e['node_from'], e['node_to']) in edge_dict:
            raise AssertionError("Multiple edges from {} to {}".format(
                e['node_from'], e['node_to']))
        if e['value'] == 0:
            raise AssertionError("Edge from {} to {} has value 0".format(
                e['node_from'], e['node_to']))

        edge_dict[(e['node_from'], e['node_to'])] = e['value']
    return edge_dict


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


def get_network_iterator():
    page = 0
    network_page = mangal_request("network", "page", page)
    while network_page:
        for net in network_page:
            yield net
        page += 1
        network_page = mangal_request("network", "page", page)


if __name__ == "__main__":
    acc = 0
    for net in get_network_iterator():
        if "poll" in net['description']:
            construct_network_excel(net['id'])
            acc += 1
    print(acc)
