import random
import xml.etree.ElementTree as elmt
from concurrent.futures import ThreadPoolExecutor

import requests as req

"""
    As it turns out, mangal already did some heavy lifting for us, there are 'taxonomy'
    entries in the database, linking to one of the following taxonomy services:
    ITIS, BOLD, EOL, NCBI, COL, GBIF. Each service has a varying percentage of coverage
    for different areas.
    
    As such, we'll implement taxonomy requests for all services both by id and my species/genus names, 
    and run them when necessary from a thread pool.
"""

DELIMETER = {'u': "_", 'p': '+'}
ACCEPTED_KINGDOMS = {'Animalia', 'Plantae'}  # added to disallow special cases like Chromista
NCBI_REQUEST_THRESHOLD = 5

MANGAL_TAXONOMY_URL = "https://mangal.io/api/v2/taxonomy/"
TAXONOMY_DB = {
    'tsn': "https://www.itis.gov/ITISWebService/jsonservice/",
    'bold': "https://v3.boldsystems.org/index.php/API_Tax/",
    'ncbi': "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    'col': "http://webservice.catalogueoflife.org/col/",
    'gbif': "https://api.gbif.org/v1/species/"
}


def taxonomy_request(request):
    """
    Base function for all web service requests. Tries to get some info from a given URL,
    raises a connection error if response is a faulty connection.
    """
    response = req.get(request)
    if not response.ok:
        raise ConnectionError("request '{0}' failed with error code {1}.".format(request, response.status_code))
    return response


def request_by_id(db, species_id):
    """
    Switch function to getting taxonomy data by id. Given a species id, and a database
    this id originates from, adds the appropriate db specific request string, and returns it's
    response (xml for NCBI, json for all else).
    :exception EoL is inconsistent in it's ids, and calling it will raise an implementation
    error
    """
    if db == 'ncbi':  # special case NCBI XML response
        id_specifier = "efetch.fcgi?db=taxonomy&id={0}".format(species_id)
        return taxonomy_request(TAXONOMY_DB['ncbi'] + str(id_specifier))
    if db == 'eol':
        # special case EOL response inconsistent with mangal taxonomy IDs
        # more info https://eol.org/docs/what-is-eol/classic-apis
        raise NotImplementedError

    if db == 'tsn':  # tsn is the name of ITIS ids
        id_specifier = "getKingdomNameFromTSN?tsn={0}".format(species_id)
    if db == 'bold':
        id_specifier = "TaxonData?dataTypes=basic&taxId={0}".format(species_id)
    if db == 'col':
        id_specifier = "webservice?format=json&response=full&id={0}".format(species_id)
    if db == 'gbif':
        id_specifier = species_id
    return taxonomy_request(TAXONOMY_DB[db] + str(id_specifier)).json()


def request_by_name(db, *args):
    """
    Switch function to getting taxonomy data by name. Given a species name, and a database
    this id originates from, adds the appropriate db specific request string, and returns it's
    response.
    :exception EoL is inconsistent in it's ids, and calling it will raise an implementation
    error;. Furthermore, ITIS and NCBI require an additional web service request to get to
    the actual taxonomy info, and deemed not worth implementing (will also raise an
    implementation error)
    """
    if db == 'tsn':
        raise NotImplementedError
    if db == 'eol':
        raise NotImplementedError
    if db == 'ncbi':
        raise NotImplementedError

    name_tokens = []
    for arg in args:
        for t in arg.split(" "):
            name_tokens.append(t)
    if db == 'bold':
        name = DELIMETER['p'].join(name_tokens)
        match_specifier = "TaxonSearch?fuzzy=false&taxName={0}".format(name)
    if db == 'col':
        name = DELIMETER['p'].join(name_tokens)
        match_specifier = "webservice?format=json&response=full&name={0}".format(name)
    if db == 'gbif':
        name = DELIMETER['u'].join(name_tokens)
        match_specifier = "match?rank=species&strict=false&name={0}".format(name)
    return taxonomy_request(TAXONOMY_DB[db] + match_specifier).json()


def itis_parse_response(response_json):
    return response_json["kingdomName"]


def bold_parse_response(response_json):
    if not response_json:
        return ""

    if "taxid" in response_json.keys():
        try:
            return response_json["tax_division"]
        except KeyError:
            return ""

    kingdom_candidates = dict()
    for _, taxon in response_json.items():
        try:
            taxon_kingdom = taxon["tax_division"]
            kingdom_candidates[taxon_kingdom] = kingdom_candidates.get(taxon_kingdom, 0) + 1
        except (KeyError, TypeError):
            continue

    if len(kingdom_candidates) > 1:
        print("BOLD: multiple possible kingdoms for species: {0}. returning none".format(kingdom_candidates))
        return ""

    selected_kingdom = ("", 0)
    for kingdom, votes in kingdom_candidates.items():
        selected_kingdom = selected_kingdom if selected_kingdom[1] > votes else (kingdom, votes)
    return selected_kingdom[0]


def eol_parse_response(response):
    raise NotImplementedError


def parse_kingdom_from_lineage(lineage):
    # NCBI uses a different nomenclature then other DBs, opting for
    # Viridoplantae (green plants) instead of plantae and
    # Metazoa instead of Animalia
    if "Metazoa" in lineage or "Animalia" in lineage:
        return "Animalia"
    elif "Viridiplantae" in lineage or "Plantae" in lineage:
        return "Plantae"
    else:
        return ""


def ncbi_parse_response(response_xml):
    root = elmt.fromstring(response_xml.text)
    taxon = root[0]  # currently ncbi requests are only by id, so only one taxon returned

    for child in taxon:
        if child.tag == "Lineage":
            return parse_kingdom_from_lineage(child.text)
    return ""


def col_parse_response(response_json):
    if response_json["number_of_results_returned"] <= 1:
        try:
            return response_json["results"][0]["classification"][0]["name"]
        except KeyError:
            return ""

    # multiple possible results in name search
    kingdom_candidates = dict()
    for result in response_json["results"]:
        try:
            result_kingdom = result["classification"][0]["name"]
            kingdom_candidates[result_kingdom] = kingdom_candidates.get(result_kingdom, 0) + 1
        except KeyError:
            pass

    if len(kingdom_candidates) > 1:
        print("CoL: multiple possible kingdoms for species: {0}. returning none".format(kingdom_candidates))
        return ""

    selected_kingdom = ("", 0)
    for kingdom, votes in kingdom_candidates.items():
        selected_kingdom = selected_kingdom if selected_kingdom[1] > votes else (kingdom, votes)
    return selected_kingdom[0]


def gbif_parse_response(response_json):
    try:
        return response_json["kingdom"]
    except KeyError:
        return ""


def parse_response(db, response):
    parsers = {
        'tsn': itis_parse_response,
        'bold': bold_parse_response,
        'ncbi': ncbi_parse_response,
        'col': col_parse_response,
        'gbif': gbif_parse_response
    }
    return parsers[db](response)


def get_mangal_taxonomy_data(tax_info):
    if type(tax_info) is int:
        tax_info = taxonomy_request(MANGAL_TAXONOMY_URL + tax_info).json()
    id_dict = dict()
    for db in TAXONOMY_DB:
        if tax_info[db] is not None:
            id_dict[db] = tax_info[db]
    return id_dict


def choose_db(db_counts, possible_dbs=None):
    possible_dbs = list(possible_dbs.keys() if type(possible_dbs) is dict else db_counts.keys())
    while True:
        chosen_db = random.choice(possible_dbs)
        if chosen_db != 'ncbi' or db_counts[chosen_db] < NCBI_REQUEST_THRESHOLD:
            break

    db_counts[chosen_db] = db_counts[chosen_db] + 1
    return chosen_db


def get_node_kingdom(tax_info=None, node_id=-1, node_kingdoms={}, db_counts=[]):
    if type(tax_info) is not str:
        tax_info = get_mangal_taxonomy_data(tax_info)

    db = choose_db(db_counts, tax_info)
    try:
        if type(tax_info) is dict:
            kingdom = parse_response(db, request_by_id(db, tax_info[db]))
        else:
            kingdom = parse_response(db, request_by_name(db, tax_info))
    except (ConnectionError, NotImplementedError) as e:
        db_counts[db] = db_counts[db] - 1
        node_kingdoms[node_id] = ""
        return

    db_counts[db] = db_counts[db] - 1
    if kingdom not in ACCEPTED_KINGDOMS:
        node_kingdoms[node_id] = ""
    node_kingdoms[node_id] = kingdom


def get_nodelist_kingdoms(nodes):
    mangal_taxonomy_info = dict()
    for node in nodes:
        if 'taxonomy' in node.keys() and node['taxonomy'] is not None:
            mangal_taxonomy_info[node['id']] = node['taxonomy']
        elif 'taxonomy_id' in node.keys() and node['taxonomy_id'] is not None:
            mangal_taxonomy_info[node['id']] = node['taxonomy_id']
        else:
            mangal_taxonomy_info[node['id']] = node['original_name']

    node_kingdoms = dict()
    db_concurrent_users = {db: 0 for db in TAXONOMY_DB.keys()}
    with ThreadPoolExecutor(max_workers=5) as executor:
        for tid, tax_info in mangal_taxonomy_info.items():
            executor.submit(get_node_kingdom, tax_info=tax_info, node_id=tid,
                            node_kingdoms=node_kingdoms, db_counts=db_concurrent_users)
    print(db_concurrent_users)  # make sure all web service requests were resolved
    return node_kingdoms


if __name__ == "__main__":
    kingdom_per_db = dict()

    for db in TAXONOMY_DB:
        try:
            kingdom_per_db[db] = parse_response(db, request_by_name(db, "Hylaeus stevensi"))
        except NotImplementedError:
            kingdom_per_db[db] = "NotImplemented"
        except (KeyError, AttributeError):
            kingdom_per_db[db] = "ParsingError"
        except ConnectionError as error:
            raise error("Could not connect to database {0}".format(db))

    print(kingdom_per_db)
