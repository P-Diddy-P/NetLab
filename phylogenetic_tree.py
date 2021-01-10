from sys import argv
import pathlib
from os import mkdir

from Bio import Phylo
from Bio.Phylo import BaseTree

from polyploid_species import PolyploidDictionary
from duplicate_finder import species_pattern


def get_phylo_options(tax_name, with_genus=False):
    def camel(s):
        return s[0].upper() + s[1:].lower()
    spec_pat = species_pattern(tax_name)
    spec_options = []

    gl, gc = spec_pat.genus.lower(), camel(spec_pat.genus)
    if len(spec_pat.species) > 0:
        sl, sc = spec_pat.species.lower(), camel(spec_pat.species)
        spec_options.extend(
            ['_'.join([gc, sl]), '_'.join([gl, sl])]
        )

    if with_genus:  # if no species is found, try searching by genus only
        spec_options.extend([gc, gl])
    return spec_options


def collapse_paths(paths):
    """
    Collapses the first node of each path as long as the first node of all paths is the same,
    and each path has more than one node in it.
    """
    if any([len(path) < 2 for path in paths]) or not len(paths):
        return

    sample_next = paths[0][0].name, paths[0][0].branch_length
    add_length = 0
    while all((path[0].name, path[0].branch_length) == sample_next for path in paths):
        add_length += sample_next[1]
        for path in paths:
            path.pop(0)
        sample_next = paths[0][0].name, paths[0][0].branch_length

        if any([len(path) < 2 for path in paths]) or not len(paths):
            break

    for path in paths:
        path[0].branch_length += add_length


def reconstruct_tree_base(root_name, tree_paths, outpath, missing):
    """
    Recursively build a new tree containing only collapsed ;param tree_paths;
    starting from ;param root_name;. The result is written to ;param outpath;
    with an auxillary file of ;param missing; species.
    """
    def reconstruct_tree_rec(root_clade, clade_paths):
        if not clade_paths:
            return

        collapse_paths(clade_paths)
        paths_per_clade = dict()
        for path in clade_paths:
            path_start = path.pop(0)
            next_node = (path_start.name, path_start.branch_length)
            next_node_paths = paths_per_clade.get(next_node, [])
            if len(path):
                next_node_paths.append(path)
            if next_node not in paths_per_clade:
                paths_per_clade[next_node] = next_node_paths

        for child_data, child_paths in paths_per_clade.items():
            child_name, child_length = child_data
            child_clade = BaseTree.Clade(name=child_name, branch_length=child_length)
            root_clade.clades.append(child_clade)
            reconstruct_tree_rec(child_clade, child_paths)

    root = BaseTree.Clade(name=root_name)
    reconstruct_tree_rec(root, tree_paths)
    Phylo.write(root, outpath.joinpath('tree'), 'newick')
    with open(outpath.joinpath('missing'), 'w') as fd:
        fd.write(str(missing))


def generate_sparse_tree(tree, species_of_interest, outpath):
    """
    Creates a restriction of ;param tree; so that it will contain only
    bifurcating nodes between ;param species_of_interest; and their LCA.
    The resulting tree is written to ;param outpath; with an auxillary file
    of species not found in ;param tree; from ;param species_of_interest;.
    """
    # First, get a list of clades to operate on when finding lca
    tree_clades = dict()
    missing_clades = set()
    for spec in species_of_interest:
        # use only major_name, the other one found to be unused on a sample of 500 species.
        # then get the first clade object (assume only one clade is returned for a species name).
        major_name, minor_name = get_phylo_options(spec)
        try:
            major_clade = [e for e in tree.find_clades(major_name)][0]
            tree_clades[spec] = major_clade
        except IndexError:
            missing_clades.add(spec)

    common_ancestor = tree.common_ancestor(tree_clades.values())
    clade_paths = dict()
    for spec, clade in tree_clades.items():
        clade_paths[spec] = common_ancestor.get_path(clade)

    try:
        mkdir(outpath)
    except FileExistsError:
        pass  # don't care about previously created directories, they will be re-written
    reconstruct_tree_base(common_ancestor.name, list(clade_paths.values()), outpath, missing_clades)


if __name__ == "__main__":
    tree_path = pathlib.Path(argv[1]).joinpath('GBMB.tre')
    polyploid_path = pathlib.Path(argv[2])

    polydict = PolyploidDictionary(polyploid_path)
    gbtree = Phylo.read(tree_path, 'newick')
