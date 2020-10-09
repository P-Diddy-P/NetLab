from sys import argv
import pathlib
import os
import re

import pandas as pd

SUSPECT_DUPLICATES = ['M_PL_001', 'M_PL_002', 'M_PL_003',
                      'arroyo_1982_19810301_949', 'arroyo_1982_19810301_950',
                      'arroyo_1982_19810301_951']


def iter_sources(sources):
    for src in sources:
        for file in os.listdir(src):
            if file.endswith('csv'):
                yield src.joinpath(file)


def extract_networks(source_directories):
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


def similar_pattern(s1, s2):
    raise NotImplementedError


def get_similar_nodes(nodes1, nodes2):
    similar_pairs = set()
    tagged1, tagged2 = set(), set()

    for n1 in nodes1:
        if n1 in tagged1:
            continue

        for n2 in nodes2:
            if n2 in tagged2:
                continue

            if similar_pattern(n1, n2):
                similar_pairs.add((n1, n2))
                tagged1.add(n1)
                tagged2.add(n2)
                break

    return similar_pairs


def compare_nodes(nodes1, nodes2, fraction=0.05):
    mean_measure = ((len(nodes1) + len(nodes2)) / 2) * fraction
    print(f'mean measure of both networks by fraction {fraction} is {mean_measure}')
    print(f'node diff is {abs(len(nodes1) - len(nodes2))}, and '
          f'{abs(len(nodes1) - len(nodes2)) / mean_measure} in mean measure')

    similar = set(nodes1).intersection(set(nodes2))
    print(f'{len(similar)} similar nodes:\n{similar}')
    print(f'{len(set(nodes1) - similar)} non-similars in net1:\n{sorted(list(set(nodes1) - similar))}')
    print(f'{len(set(nodes2) - similar)} non-similars in net1:\n{sorted(list(set(nodes2) - similar))}')


def special_pattern(string):
    special_chars = r'[^a-zA-Z]'
    for word in string.split(' '):
        if re.search(special_chars, word):
            return True
    return False


if __name__ == "__main__":
    # TODO currently doesn't parse IWDB networks, as they have no set format.
    # TODO add IWDB (argv[1]) once it's directory is in order.
    # networks = extract_networks([pathlib.Path(s) for s in argv[2:]])
    sources = [pathlib.Path(argv[2])]

    for net_name, network in extract_networks(sources).items():
        plant_sp = 0
        pol_sp = 0
        print(f"============================================\n"
              f"searching {net_name}\n")
        for plant_name in network.index:
            if special_pattern(plant_name):
                plant_sp += 1
                print(plant_name)

        for pol_name in network.columns:
            if special_pattern(pol_name):
                pol_sp += 1
                print(pol_name)
        print(f"\ntotal {plant_sp} special plants || {pol_sp} special pollinators")
        print("============================================\n")
