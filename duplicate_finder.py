from sys import argv
from collections import namedtuple
import pathlib
import os
import re

import pandas as pd

SUSPECT_DUPLICATES = ['M_PL_001', 'M_PL_002', 'M_PL_003',
                      'arroyo_1982_19810301_949', 'arroyo_1982_19810301_950',
                      'arroyo_1982_19810301_951']

SPECIES_SPATTERN = r' (?P<base>sp\.{0,2}|n\.i\.) ?' \
                   r'(?P<spindex>[0-9]*) ?' \
                   r'(?P<spnet>M_PL_[0-9]{3}(_[0-9]{0,3})?)? ?' \
                   r'(?P<spextra>.*)'

SpeciesIdentifier = namedtuple('SpeciesIdentifier',
                               ['genus', 'index', 'network', 'extra', 'origin'])


def iter_sources(sources):  # TODO move to mainframe when done!
    for src in sources:
        for file in os.listdir(src):
            if file.endswith('csv'):
                yield src.joinpath(file)


def extract_networks(source_directories):   # TODO move to mainframe when done!
    all_networks = dict()
    duplicate_names = []
    for net_file in iter_sources(source_directories):
        if net_file.stem not in all_networks:
            all_networks[net_file.stem] = pd.read_csv(net_file, index_col=0)
        else:
            duplicate_names.append(net_file.stem)

    if duplicate_names:
        print(f"found multiple networks with names: {', '.join(duplicate_names)}")
    return all_networks


def species_pattern(species):
    """
    Parses a names of unidentified species to the following tokens:
    - genus
    - index (for when there are several unidentified species)
    - network identifier (special for WoL networks with M_PL_[XXX])
    - Extra information (usually a few capitalized letters)
    This breaking down of unknown species will allow for comparison by subsets
    or specific criteria.
    :param species: string with the species name as given in the network.
    :return: a regularized namedtuple of the species name
    """
    normal_species = SpeciesIdentifier(
        genus='', index='', network='', extra='', origin=species
    )
    if len(species.split(' ')) == 2 and not \
            any([re.search(r'(^sp$|[^a-zA-Z])', w) for w in species.split(' ')]):
        # Don't consider species with "valid" structure of `genus species`
        return normal_species

    match = re.search(SPECIES_SPATTERN, species)
    if not match:
        return normal_species

    return SpeciesIdentifier(
        genus=species.split(' ')[0],
        index=match.group('spindex'),
        network='' if not match.group('spnet') else match.group('spnet'),
        extra=match.group('spextra'),
        origin=species
    )


def parameterized_compare(a, b, *args):
    """
    Compares two elements: a and b, according to a list
    of parameters, given as args.
    :param a: first comparand.
    :param b: second comparand.
    :param args: parameters to compare by. If no parameters
    are given, compare a == b.
    :return: True if a and b have the same parameters.
    """
    if not args:
        return a == b
    return all([getattr(a, arg) == getattr(b, arg) for arg in args])


def map_nodeset_by_parameters(set1, set2, map12, map21, *args):
    for cmpa in set1:
        if cmpa in map12:
            pass

        for cmpb in set2:
            if parameterized_compare(cmpa, cmpb, *args) \
                    and not map21.get(cmpb, None):
                map12[cmpa.origin] = cmpb.origin
                map21[cmpb.origin] = cmpa.origin


def map_nodeset(set1, set2):
    """
    Given two sets of network nodes labeled by species names, maps
    the names between the nodes in each set, first by full text comparison,
    then by a diminishing number of parameters using the species pattern.
    Assumes all node names in a network are unique.
    :return: quadruple of the following:
    1. mapping from set1 to set2 nodes.
    2. reverse mapping from set2 to set1 nodes.
    3. nodes without mapping in set1.
    4. nodes without mapping in set2.
    Note that nodes without mapping are mapped to None in 1. and 2.
    """
    map1to2 = dict()
    map2to1 = dict()

    pset1 = [species_pattern(sp_string) for sp_string in set1]
    pset2 = [species_pattern(sp_string) for sp_string in set2]
    map_nodeset_by_parameters(pset1, pset2, map1to2, map2to1, 'origin')
    print(map1to2)
    print(map2to1)


def compare_networks(net1, net2):
    """
    Takes two networks and compares them according to matching in
    plant/pollinator/interaction elements and set sizes.
    :return: k-tuple containing two mappings from net1 vertices
     to net2 and vice versa, as well as several similarity/coverage
     parameters. If no adequate mappings are found, they are replaced
     by None values.
    """
    raise NotImplementedError


def count_unidentified(net_table):
    """
    Counts the number of COMPLETELY unidentified species of plants and pollinators
    in a network. Assumes such species are named `unidentified` in the network,
    and are wholly unknown (i.e. the genus is also unknown).
    :param net_table: pandas table of a pollinator network
    :return: tuple of unknown plant and pollinator species.
    """
    return sum(1 if p.lower().startswith('unidentified') else 0 for p in net_table.index), \
        sum(1 if p.lower().startswith('unidentified') else 0 for p in net_table.columns)


if __name__ == "__main__":
    # Decided method for duplicate finding:
    # remove unknown species from both plant and pollinator lists
    # (drop networks with over X% unknown species of plants or pollinators).
    # After removing (ENTIRELY) unknown species, preprocess both plant
    # and pollinator list by formalizing names to some format of required
    # parameters for comparison, and non required.
    # Then, set intersect plant and pollinator lists for networks, and compare
    # all other parts manually by format (consider comparing all by format).
    # If all goes well with nodelist comparisons, compare interactions
    # (for efficient interaction comparison, keep base species strings to access
    # pd cells easily).

    sources = [pathlib.Path(arg) for arg in argv[2:]]
    net_pd = extract_networks(sources)
    map_nodeset(net_pd[SUSPECT_DUPLICATES[2]].index, net_pd[SUSPECT_DUPLICATES[5]].index)

    """
    for net_name, network in extract_networks(sources).items():
        plant_sp = 0
        pol_sp = 0
        print(f"============================================\n"
              f"searching {net_name}\n")
        for plant_name in network.index:
            if pattern := species_pattern(plant_name):
                plant_sp += 1
                print(pattern)

        for pol_name in network.columns:
            if pattern := species_pattern(pol_name):
                pol_sp += 1
                print(pattern)
        print(f"\ntotal {plant_sp} special plants || {pol_sp} special pollinators")
        print(count_unidentified(network))
    """
