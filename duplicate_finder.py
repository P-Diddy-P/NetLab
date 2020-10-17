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
                               ['genus', 'species', 'index', 'network', 'extra', 'origin'])


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


def name_missing_structure(tax_name):
    """
    Checks whether the taxon name follows the pattern expected
    from a taxon with missing species, that is: genus (might be shorthand),
    followed by an unknown species identifier.
    :param tax_name: taxon name to check.
    :return: True if the taxon name follows an unknown species pattern,
    otherwise False.
    """
    query_word = tax_name.split(' ')[1]
    return re.match(r'(^sp$|n\.i\.|^sp[^a-z])', query_word)


def species_pattern(tax_name):
    """
    Parses a names of unidentified species to the following tokens:
    - genus
    - index (for when there are several unidentified species)
    - network identifier (special for WoL networks with M_PL_[XXX])
    - Extra information (usually a few capitalized letters)
    This breaking down of unknown species will allow for comparison by subsets
    or specific criteria.
    :param tax_name: string with the species name as given in the network.
    :return: a regularized namedtuple of the species name
    """
    tax_words = tax_name.split(' ')
    normal_structure = SpeciesIdentifier(
        genus=tax_words[0],
        species=tax_words[1],
        index='', network='', origin=tax_name,
        extra=tax_words[2:] if len(tax_words) > 2 else ''
    )
    if not name_missing_structure(tax_name):
        return normal_structure

    match = re.search(SPECIES_SPATTERN, tax_name)
    if not match:
        return normal_structure

    return SpeciesIdentifier(
        genus=tax_name.split(' ')[0],
        species='',
        index=match.group('spindex'),
        network='' if not match.group('spnet') else match.group('spnet'),
        extra=match.group('spextra'),
        origin=tax_name
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
        if cmpa.origin in map12:
            continue

        for cmpb in set2:
            if cmpb.origin in map21:
                continue

            if parameterized_compare(cmpa, cmpb, *args) \
                    and not map21.get(cmpb, None):
                # print("mapping ", cmpa.origin, " to ", cmpb.origin, " by ", args)
                map12[cmpa.origin] = cmpb.origin
                map21[cmpb.origin] = cmpa.origin
                break


def map_non_matches(key_set, mapping):
    """
    Completes a mapping to contain all keys in key_set, with missing
    keys mapped to None.
    :param key_set: a set of elements which should be keys
    in mapping.
    :param mapping: dictionary mapping to complete.
    :return: None, but mutates mapping patameter.
    """
    for key in key_set:
        if key not in mapping.keys():
            mapping[key] = None


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
    nonspecific_parameters = ['genus', 'index', 'network', 'extra']
    map1to2, map2to1 = dict(), dict()

    pset1 = [species_pattern(sp_string) for sp_string in set1]
    pset2 = [species_pattern(sp_string) for sp_string in set2]
    map_nodeset_by_parameters(pset1, pset2, map1to2, map2to1, 'origin')
    for stop in range(len(nonspecific_parameters), 1, -1):
        map_nodeset_by_parameters(pset1, pset2, map1to2, map2to1,
                                  *nonspecific_parameters[:stop])
    map_nodeset_by_parameters(pset1, pset2, map1to2, map2to1, 'species')

    unmapped1, unmapped2 = set1 - set(map1to2.keys()), set2 - set(map2to1.keys())
    map_non_matches(set1, map1to2)
    map_non_matches(set2, map2to1)
    return map1to2, map2to1, unmapped1, unmapped2


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


def compare_networks(net1, net2):
    """
    Takes two networks and compares them according to matching in
    plant/pollinator/interaction elements and set sizes.
    :return: k-tuple containing two mappings from net1 vertices
     to net2 and vice versa, as well as several similarity/coverage
     parameters.
    """
    raise NotImplementedError


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
    map_nodeset(
        set(net_pd[SUSPECT_DUPLICATES[2]].index),
        set(net_pd[SUSPECT_DUPLICATES[5]].index)
    )

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
