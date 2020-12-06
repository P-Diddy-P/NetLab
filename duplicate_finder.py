from sys import argv
from collections import namedtuple
import random
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


def name_missing_structure(tax_name):
    """
    Checks whether the taxon name follows the pattern expected
    from a taxon with missing species, that is: genus (might be shorthand),
    followed by an unknown species identifier.
    :param tax_name: taxon name to check.
    :return: True if the taxon name follows an unknown species pattern,
    otherwise False.
    """
    try:
        query_word = tax_name.split(' ')[1]
        return re.match(r'(^sp$|n\.i\.|^sp[^a-z])', query_word)
    except IndexError:  # some tax names are just one word, can't do much with that
        return None


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
    try:
        tax_words = tax_name.split(' ')
        normal_structure = SpeciesIdentifier(
            genus=tax_words[0],
            species=tax_words[1] if len(tax_words) > 1 else '',
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
    except IndexError:
        print(f"CAN'T FIND SPECIES PATTERN FOR: {tax_name}")
        raise


def parameterized_compare(a, b, allow_empty, *args):
    """
    Compares two elements: a and b, according to a list
    of parameters, given as args.
    :param a: first comparand.
    :param b: second comparand.
    :param allow_empty: allows comparing empty parameters. If False,
    then (a.arg == b.arg && a.arg == '') will return False
    :param args: parameters to compare by. If no parameters
    are given, compare a == b.
    :return: True if a and b have the same parameters.
    """
    if not args:
        return a == b
    return all([getattr(a, arg) == getattr(b, arg) and
                (getattr(a, arg) != '' or allow_empty)
                for arg in args])


def map_nodeset_by_parameters(set1, set2, map12, map21, cmp_empty, *args):
    """
    Takes unmapped elements from set1 and set2 and maps between them by
    the given comparison parameters.
    :param set1: first set of mapping elements.
    :param set2: second set of mapping elements.
    :param map12: mapping from set1 elements to set2.
    :param map21: mapping from set2 elements to set1.
    :param cmp_empty: consider empty parameters comparable.
    :param args: comparison parameters.
    :return: None, but mutates map parameters.
    """
    for cmpa in set1:
        if cmpa.origin in map12:
            continue

        for cmpb in set2:
            if cmpb.origin in map21:
                continue

            if parameterized_compare(cmpa, cmpb, not cmp_empty, *args) \
                    and not map21.get(cmpb.origin, None):
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
    map_nodeset_by_parameters(pset1, pset2, map1to2, map2to1, True, 'origin')
    for stop in range(len(nonspecific_parameters), 1, -1):
        map_nodeset_by_parameters(pset1, pset2, map1to2, map2to1, False,
                                  *nonspecific_parameters[:stop])
    map_nodeset_by_parameters(pset1, pset2, map1to2, map2to1, True, 'species')

    unmapped1, unmapped2 = set1 - set(map1to2.keys()), set2 - set(map2to1.keys())
    map_non_matches(set1, map1to2)
    map_non_matches(set2, map2to1)
    return map1to2, map2to1, unmapped1, unmapped2


def compare_network_interactions(net1, net2, plant_map12, plant_map21,
                                 pol_map12, pol_map21, compare_value=False,
                                 count_nomaps=False):
    """
    Compares the ineraction matrices between two networks using the mappings
    between both networks, and counts the number of discrepancies (i.e. interaction
    between two species in one network, but not in the other.
    :param count_nomaps: Count missing interaction when one of the nodes has no
    mapping to the other network.
    :param compare_value: Count interactions with different values as "missing".
    :return: A quadruple containing:
        1. Number of missing interactions in net1.
        2. Number of missing interactions in net2.
        3. Total number of interactions in net1.
        4. Total number of interactions in net2.
    """
    # go over both networks, when going over network1, count interactions to
    # net1 total and missing interactions to net2 missing. Likewise for net2.
    net1_total, net2_total = 0, 0
    net1_missing, net2_missing = 0, 0

    for plant_name, plant_row in net1.iterrows():
        for pol_name, pol_interaction in plant_row.iteritems():
            if pol_interaction == 0:
                continue

            net1_total += 1
            mapped_plant = plant_map12[plant_name]
            mapped_pol = pol_map12[pol_name]
            # print(f"({plant_name},{pol_name}) => ({mapped_plant}, {mapped_pol})")

            if not (mapped_pol and mapped_plant):
                net2_missing += (1 if count_nomaps else 0)
                continue

            if (net2.loc[mapped_plant, mapped_pol] != pol_interaction and compare_value) or \
                    not net2.loc[mapped_plant, mapped_pol]:
                net2_missing += 1

    for plant_name, plant_row in net2.iterrows():
        for pol_name, pol_interaction in plant_row.iteritems():
            if pol_interaction == 0:
                continue

            net2_total += 1
            mapped_plant = plant_map21[plant_name]
            mapped_pol = pol_map21[pol_name]
            if not (mapped_pol and mapped_plant):
                net1_missing += (1 if count_nomaps else 0)
                continue

            if (net1.loc[mapped_plant, mapped_pol] != pol_interaction and compare_value) or \
                    not net1.loc[mapped_plant, mapped_pol]:
                net1_missing += 1

    return net1_missing, net2_missing, net1_total, net2_total


def compare_networks(net1, net2, ct=0.05, log=False):
    """
    Takes two networks and compares them according to matching in
    plant/pollinator/interaction elements and set sizes.
    :param ct: comparison threshold between networks (multiplied by
    mean plant/pollinator/interaction number).
    :param log: print the whether the networks are duplicate,
    and why not.
    :return: True if networks are duplicates, False otherwise.
    """
    plantsize1, polsize1 = len(net1.index), len(net1.columns)
    plantsize2, polsize2 = len(net2.index), len(net2.columns)
    plant_threshold = (plantsize1 + plantsize2) / 2 * ct
    pol_threshold = (polsize1 + polsize2) / 2 * ct

    if abs(plantsize1 - plantsize2) > plant_threshold:
        if log:
            print(f"-->Networks are not duplicates due to different plant sizes: "
                  f"{plantsize1} || {plantsize2}")
        return False
    if abs(polsize1 - polsize2) > pol_threshold:
        if log:
            print(f"-->Networks are not duplicate due to different pollinator sizes: "
                  f"{polsize1} || {polsize2}")
        return False

    plants1to2, plants2to1, plants1_unmapped, plants2_unmapped = map_nodeset(
        set(net1.index), set(net2.index))
    pols1to2, pols2to1, pols1_unmapped, pols2_unmapped = map_nodeset(
        set(net1.columns), set(net2.columns))

    if len(plants1_unmapped | plants2_unmapped) > plant_threshold:
        if log:
            print(f"-->Networks are not duplicate due to plant mapping failure"
                  f" with unmapped {len(plants1_unmapped)} || {len(plants2_unmapped)}")
        return False
    if len(pols1_unmapped | pols2_unmapped) > pol_threshold:
        if log:
            print(f"-->Networks are not duplicate due to pollinator mapping failure"
                  f" with unmapped {len(pols1_unmapped)} || {len(pols2_unmapped)}")
        return False

    net1_missing, net2_missing, net1_total, net2_total = compare_network_interactions(
        net1, net2, plants1to2, plants2to1, pols1to2, pols2to1)
    interaction_threshold = (net1_total + net2_total) / 2 * ct
    if (net1_missing + net2_missing) > interaction_threshold:
        if log:
            print(f"-->Networks are not duplicate due to different interactions with missing "
                  f"{net1_missing} || {net2_missing}")
        return False

    if log:
        print("-->duplicate networks!")
    return True


def drop_network(net1, net2):
    """
    Compares networks by plantsize, then by polsize, then randomly
    :return: True if should drop net2, False if should drop net1
    """
    plantsize1, polsize1 = len(net1.index), len(net1.columns)
    plantsize2, polsize2 = len(net2.index), len(net2.columns)

    if plantsize1 != plantsize2:
        return plantsize1 > plantsize2

    if polsize1 != polsize2:
        return polsize1 > polsize2

    return random.random() > 0.5


def compare_all_networks(networks, threshold, drop_early=frozenset(), log=False):
    net_keys = list(networks.keys())
    duplicate_networks = set(drop_early)

    for i in range(len(net_keys)):
        net_i = net_keys[i]
        if net_i in duplicate_networks:
            break

        for j in range(i + 1, len(net_keys)):
            net_j = net_keys[j]
            if net_j in duplicate_networks:
                break

            if log:
                print(f"\ncomparing: [{net_i}] || [{net_j}]")

            if compare_networks(networks[net_i], networks[net_j],
                                ct=threshold, log=log):
                drop_j = drop_network(networks[net_i], networks[net_j])
                duplicate_networks.add(net_j if drop_j else net_i)

                if not drop_j:
                    break

    return duplicate_networks


if __name__ == "__main__":
    print("removed network utility function to mainframe.")
    """
    sources = [pathlib.Path(arg) for arg in argv[2:]]
    net_pd = extract_networks(sources)
    unknown_networks = clean_networks(net_pd)

    net_keys = list(net_pd.keys())
    for i in range(len(net_keys)):
        for j in range(i + 1, len(net_keys)):
            print(f"\ncomparing: [{net_keys[i]}] || [{net_keys[j]}]")
            compare_networks(net_pd[net_keys[i]], net_pd[net_keys[j]], ct=0.05)
    """
