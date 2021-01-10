import unittest
import io

from Bio import Phylo

from phylogenetic_tree import *


class TreeReconstruction(unittest.TestCase):
    full_tree = "((((B:1.2,C:1.7)A:2.34)IA2:0.44)IA1:0.12,(((((L:0.74)IL3:0.07)IL2:0.03)IL1:0.14,((F:0.23)IF1:1.01,(G:0.4)IG1:0.98)E:0.52)D:1.08)ID1:0.79)R:0.0;"
    reduced_tree = "((B:1.20000,C:1.70000)A:2.90000,(L:0.98000,(F:1.24000,G:1.38000)E:0.52000)D:1.87000)R:0.00000;\n"

    def test_tree_reduction(self):
        full_tree = Phylo.read(io.StringIO(self.full_tree), 'newick')

        leaves = full_tree.get_terminals()
        lca = full_tree.common_ancestor(leaves)
        leaf_paths = [lca.get_path(leaf) for leaf in leaves]
        reconstructed_tree = reconstruct_tree_base(lca.name, leaf_paths)

        reconstructed_string = io.StringIO()
        Phylo.write(reconstructed_tree, reconstructed_string, 'newick')
        assert self.reduced_tree == reconstructed_string.getvalue()


if __name__ == "__main__":
    unittest.main()
