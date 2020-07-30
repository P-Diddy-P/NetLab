import unittest
from taxonomy_info import *


def get_kingdoms_per_db_by_id(mid):
    taxonomy_ids = get_mangal_taxonomy_data(mid)
    kingdom_per_db = dict()

    for db, tid in taxonomy_ids.items():
        kingdom_per_db[db] = parse_response(db, request_by_id(db, tid))
    return kingdom_per_db


def get_kingdoms_per_db_by_name(*args):
    kingdom_per_db = dict()

    for db in TAXONOMY_DB:
        try:
            kingdom_per_db[db] = parse_response(db, request_by_name(db, *args))
        except NotImplementedError:
            kingdom_per_db[db] = "NotImplemented"
        except (KeyError, AttributeError):
            kingdom_per_db[db] = "ParsingError"
        except ConnectionError as error:
            raise error("Could not connect to database {0}".format(db))
    return kingdom_per_db


class TaxonomyById(unittest.TestCase):
    def test_plant_taxonomy_by_id(self):
        kingdoms = get_kingdoms_per_db_by_id("2")
        self.assertEqual(kingdoms, {
            'tsn': 'Plantae',
            'bold': 'Plantae',
            'ncbi': 'Plantae',
            'col': 'Plantae',
            'gbif': 'Plantae'
        })

    def test_animal_taxonomy_by_id(self):
        kingdoms = get_kingdoms_per_db_by_id("5317")
        self.assertEqual(kingdoms, {
            'tsn': 'Animalia',
            'bold': 'Animalia',
            'ncbi': 'Animalia',
            'col': 'Animalia',
            'gbif': 'Animalia'
        })


class FullNetworkTaxonomy(unittest.TestCase):
    def test_full_network(self):
        # hijacking taxonomy_request to get the nodes of network 1295 (has only 5 nodes)
        net_nodes = taxonomy_request('https://mangal.io/api/v2/node?network_id=1295')
        node_kingdoms = get_nodelist_kingdom(net_nodes.json())
        self.assertEqual(node_kingdoms, {
            19907: 'Plantae',
            19908: 'Plantae',
            19909: 'Animalia',
            19910: 'Animalia',
            19911: 'Animalia'
        })

    def test_network_coverage(self):
        net_nodes = taxonomy_request('https://mangal.io/api/v2/node?network_id=27')
        self.assertEqual(len(net_nodes.json()), len(get_nodelist_kingdom(net_nodes.json())))


if __name__ == '__main__':
    unittest.main()
