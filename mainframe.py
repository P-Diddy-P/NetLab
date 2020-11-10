from sys import argv
import pathlib
import os

import pandas as pd

from polyploid_species import extract_polyploids
from duplicate_finder import compare_all_networks


def unidentified_rate(net_table):
    """
    Counts the number of COMPLETELY unidentified species of plants and pollinators
    in a network and divides by the total species number. Assumes such species are
    named `unidentified` in the network, and are wholly unknown (i.e. the genus is
    also unknown).
    :param net_table: pandas table of a pollinator network
    :return: tuple of unknown plant and pollinator rates.
    """
    unidentified_plants = sum(1 if p.lower().startswith('unidentified') else 0
                              for p in net_table.index)
    unidentified_pols = sum(1 if p.lower().startswith('unidentified') else 0
                            for p in net_table.columns)
    return unidentified_plants / len(net_table.index), unidentified_pols / len(net_table.columns)


def clean_networks(networks, plant_threshold=0.25, pol_threshold=1.0,
                   clean_plants=True, clean_pollinators=False):
    partial_networks = set()

    for net_name in networks.keys():
        ref_net = networks[net_name]

        try:
            missing_plants, missing_pols = unidentified_rate(ref_net)
        except:
            print(net_name, "\n", ref_net)
            raise
        known_cols = [not (c.lower().startswith('unidentified') and clean_plants)
                      for c in ref_net.columns]
        known_rows = [not (r.lower().startswith('unidentified') and clean_pollinators)
                      for r in ref_net.index]

        dropped_net = ref_net.loc[known_rows, known_cols]
        networks[net_name] = dropped_net

        if missing_plants > plant_threshold or missing_pols > pol_threshold:
            partial_networks.add(net_name)

    return partial_networks


def iter_sources(source_directories):
    for src in source_directories:
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

    return all_networks, duplicate_names


if __name__ == "__main__":
    polyploid_path = pathlib.Path(argv[1])
    network_paths = [pathlib.Path(path_str) for path_str in argv[2:]]

    definite_polyploids, suspect_polyploids = extract_polyploids(polyploid_path)
    network_tables, duplicate_network_names = extract_networks(network_paths)
    partial_networks = clean_networks(network_tables)

    duplicates = compare_all_networks(network_tables, 0.05, drop_early=partial_networks, duplicate_reason=False)
    print(len(network_tables) - len(duplicates))