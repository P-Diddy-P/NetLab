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
    if any([len(path) < 2 for path in paths]) or len(paths) < 1:
        return

    sample_next = paths[0][1].name, paths[0][1].branch_length
    add_length = 0
    while all((path[1].name, path[1].branch_length) == sample_next for path in paths):
        add_length += paths[0][0].branch_length
        for path in paths:
            path.pop(0)
        try:
            sample_next = paths[0][1].name, paths[0][1].branch_length
        except IndexError:
            break

    # the first node on each path should refer to the same clade
    paths[0][0].branch_length += add_length


def reconstruct_tree_base(root_name, tree_paths):
    """
    Recursively build a new tree containing only collapsed ;param tree_paths;
    starting from ;param root_name; and return the resulting tree.
    """
    def reconstruct_tree_rec(root_clade, clade_paths):
        if not clade_paths:
            return

        paths_per_clade = dict()
        for path in clade_paths:
            path_start = path[0]
            next_node = (path_start.name, path_start.branch_length)
            next_node_paths = paths_per_clade.get(next_node, [])
            next_node_paths.append(path)
            if next_node not in paths_per_clade:
                paths_per_clade[next_node] = next_node_paths

        for child_data, child_paths in paths_per_clade.items():
            collapse_paths(child_paths)
            for path in child_paths:
                child_data = path.pop(0)
            child_name, child_length = child_data.name, child_data.branch_length
            child_clade = BaseTree.Clade(name=child_name, branch_length=child_length)
            root_clade.clades.append(child_clade)

            child_paths = [p for p in child_paths if len(p)]
            reconstruct_tree_rec(child_clade, child_paths)

    root = BaseTree.Clade(name=root_name)
    reconstruct_tree_rec(root, tree_paths)
    return root


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
        try:
            major_name, minor_name = get_phylo_options(spec)
            major_clade = [e for e in tree.find_clades(major_name)][0]
            tree_clades[spec] = major_clade
        except (IndexError, ValueError) as e:
            missing_clades.add(spec)

    common_ancestor = tree.common_ancestor(tree_clades.values())
    clade_paths = dict()
    for spec, clade in tree_clades.items():
        clade_paths[spec] = common_ancestor.get_path(clade)
    reconstructed_tree = reconstruct_tree_base(common_ancestor.name, list(clade_paths.values()))

    try:
        mkdir(outpath)
    except FileExistsError:
        pass  # don't care about previously created directories, they will be re-written
    Phylo.write(reconstructed_tree, outpath.joinpath('tree.tre'), 'newick')
    with open(outpath.joinpath('missing'), 'w') as fd:
        fd.write(str(missing_clades))


def closest_leaf_distances(tree):
    """
    Creates a dictionary of a leaf mapped to the shortest distance between it and another leaf.
    """
    distances = {leaf: float("inf") for leaf in tree.get_terminals()}
    for leaf in tree.get_terminals():
        path = tree.get_path(leaf)
        path.insert(0, tree.root)
        height = path.pop().branch_length

        while height < distances[leaf] and len(path) > 0:
            current_parent = path.pop()
            search_leaves = current_parent.get_terminals()
            search_leaves.remove(leaf)

            for s_leaf in search_leaves:

                distance = height + current_parent.distance(s_leaf)
                if distance < distances[leaf]:
                    distances[leaf] = distance
                if distance < distances[s_leaf]:
                    distances[s_leaf] = distance

            height += current_parent.branch_length

    depth = max([tree.distance(leaf) for leaf in tree.get_terminals()])
    distances = {leaf: distances[leaf] / depth for leaf in tree.get_terminals()}
    return distances


if __name__ == "__main__":
    tree_path = pathlib.Path(argv[1]).joinpath('GBMB.tre')
    polyploid_path = pathlib.Path(argv[2])

    polydict = PolyploidDictionary(polyploid_path)
    gbtree = Phylo.read(tree_path, 'newick')
