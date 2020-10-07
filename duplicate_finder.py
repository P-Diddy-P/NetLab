from sys import argv
import pathlib
import os

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


if __name__ == "__main__":
    # TODO currently doesn't parse IWDB networks, as they have no set format.
    # TODO add IWDB (argv[1]) once it's directory is in order.
    networks = extract_networks([pathlib.Path(s) for s in argv[2:]])

    for i in range(3):
        wol_sus = networks[SUSPECT_DUPLICATES[i]]
        mangal_sus = networks[SUSPECT_DUPLICATES[i + 3]]
        print("comparing suspicious networks:\n"
              "===============================================")
        print("plants are kinda sus...\n"
              "-----------------------------------------------")
        compare_nodes(wol_sus.index, mangal_sus.index)
        print("pollinators are kinda sus...\n"
              "-----------------------------------------------")
        compare_nodes(wol_sus.columns, mangal_sus.columns)

    """
    net_keys = list(networks.keys())
    for i in range(len(net_keys)):
        for j in range(i+1, len(net_keys)):
            if check_for_duplicate(networks[net_keys[i]],
                                networks[net_keys[j]]):
                print(f"networks {net_keys[i]} and {net_keys[j]} are suspect duplicates")
    """
