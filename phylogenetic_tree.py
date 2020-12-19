from sys import argv
import pathlib

from Bio import Phylo

from polyploid_species import PolyploidDictionary


def verify_in_tree(species_of_interest, tree):
    leaves = set([leaf.name.lower() for leaf in tree.get_terminals()])
    missing = 0
    for species in species_of_interest:
        species = '_'.join(species.split(' ')).lower()
        if species not in leaves:
            print(f"{species} not found in tree!")
            missing += 1
    return missing


if __name__ == "__main__":
    tree_path = pathlib.Path(argv[1]).joinpath('ALLMB.tre')
    polyploid_path = pathlib.Path(argv[2])

    polydict = PolyploidDictionary(polyploid_path)
    gbtree = Phylo.read(tree_path, 'newick')

    verify_in_tree(polydict.iterate_polyploids(), gbtree)
