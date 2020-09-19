import unittest
from mangal_parser import *


class NetworkKingdomCompletion(unittest.TestCase):
    def test_invalid_networks(self):
        """
        Test mangal_parser raising alerts for non-bipartite networks encountered.
        In this case the parser should raise an assertion error.
        """
        faux_vertices = {
            '5': {'id': '5', 'name': 'v5', 'kingdom': ''},
            '4': {'id': '4', 'name': 'v4', 'kingdom': ''},
            '3': {'id': '3', 'name': 'v3', 'kingdom': ''},
            '2': {'id': '2', 'name': 'v2', 'kingdom': ''},
            '1': {'id': '1', 'name': 'v1', 'kingdom': 'Plantae'}
        }
        faux_edges = {
            ('1', '2'): 1,
            ('2', '3'): 1,
            ('3', '4'): 1,
            ('4', '5'): 1,
            ('1', '5'): 1
        }

        with self.assertRaises(AssertionError):
            complete_kingdoms_by_interactions(faux_vertices, faux_edges)

    def test_isolated_components(self):
        """
        Test how mangal_parser handles a graph with two components when
        there is no taxonomy information about any node in one component.
        In this case the parser should complete the kingdom in whatever components
        it can, leaving the unknown component untouched.
        """
        faux_vertices = {
            '1': {'id': '1', 'name': 'v1', 'kingdom': ''},
            '2': {'id': '2', 'name': 'v2', 'kingdom': ''},
            '3': {'id': '3', 'name': 'v3', 'kingdom': 'Plantae'},
            '4': {'id': '4', 'name': 'v4', 'kingdom': ''},
            '5': {'id': '5', 'name': 'v5', 'kingdom': ''},
            '6': {'id': '6', 'name': 'i6', 'kingdom': ''},
            '7': {'id': '7', 'name': 'i7', 'kingdom': ''}
        }
        faux_edges = {
            ('1', '2'): 1,
            ('2', '3'): 1,
            ('3', '4'): 1,
            ('4', '5'): 1,
            ('6', '7'): 1,
        }
        expected_result = {
            '1': 'Plantae',
            '2': 'Animalia',
            '3': 'Plantae',
            '4': 'Animalia',
            '5': 'Plantae',
            '6': '',
            '7': ''
        }

        complete_kingdoms_by_interactions(faux_vertices, faux_edges)
        self.assertEqual({k: v['kingdom'] for k, v in faux_vertices.items()},
                         expected_result)

    def test_multiple_cycle_completion(self):
        """
        Test mangal_parser handling of networks which require multiple rounds
        of inferring kingdoms by edge connections to classify all nodes by kingdom.
        In this case the network should be wholly classified after 4 completion cycles.
        """
        faux_vertices = {
            '5': {'id': '5', 'name': 'v5', 'kingdom': ''},
            '4': {'id': '4', 'name': 'v4', 'kingdom': ''},
            '3': {'id': '3', 'name': 'v3', 'kingdom': ''},
            '2': {'id': '2', 'name': 'v2', 'kingdom': ''},
            '1': {'id': '1', 'name': 'v1', 'kingdom': 'Plantae'}
        }
        faux_edges = {
            ('1', '2'): 1,
            ('2', '3'): 1,
            ('3', '4'): 1,
            ('4', '5'): 1
        }
        expected_result = {
            '1': 'Plantae',
            '2': 'Animalia',
            '3': 'Plantae',
            '4': 'Animalia',
            '5': 'Plantae'
        }

        complete_kingdoms_by_interactions(faux_vertices, faux_edges)
        self.assertEqual({k: v['kingdom'] for k, v in faux_vertices.items()},
                         expected_result)


if __name__ == '__main__':
    unittest.main()
