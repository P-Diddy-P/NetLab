from sys import argv
import pathlib
import random
import os

import pandas as pd

from polyploid_species import PolyploidDictionary
from duplicate_finder import compare_all_networks
from network_analysis import rank_graph, network_nodf


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
                   clean_plants=False, clean_pollinators=False):
    partial_networks = set()

    for net_name in networks.keys():
        ref_net = networks[net_name]

        try:
            missing_plants, missing_pols = unidentified_rate(ref_net)
        except Exception:
            print(net_name, "\n", ref_net)
            raise

        if clean_pollinators or clean_plants:
            known_cols = [not (c.lower().startswith('unidentified') and clean_pollinators)
                          for c in ref_net.columns]
            known_rows = [not (r.lower().startswith('unidentified') and clean_plants)
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
            duplicate_names.append(net_file)

    return all_networks, duplicate_names


def sample_networks(networks, no_pick=None, size=10):
    viable_networks = set(networks.keys())
    if no_pick:
        viable_networks = viable_networks - no_pick
    return random.sample(viable_networks, size)


def analyze_networks(networks, polyploids, measure='pg'):
    print(f"Insert network analysis pipe here")
    for name, table in networks.items():
        ranked_species = rank_graph(table, measure_by=measure)
        print(ranked_species)


if __name__ == "__main__":
    polyploid_path = pathlib.Path(argv[1])
    network_paths = [pathlib.Path(path_str) for path_str in argv[2:]]

    polydict = PolyploidDictionary(polyploid_path)
    network_tables, duplicate_network_names = extract_networks(network_paths)
    invalid_networks = clean_networks(network_tables)

    duplicates = compare_all_networks(network_tables, 0.05, drop_early=invalid_networks)
    print(f"{len(duplicates)} duplicates found")
    network_polyploids = dict()
    for name, table in network_tables.items():
        polies = {pn for pn in table.index if polydict.test_ploidy(pn)}
        if len(polies) > 0:
            network_polyploids[name] = polies

    networks_to_analyze = {name: table for name, table in network_tables.items() if
                           name in network_polyploids and
                           name not in duplicates and
                           name not in invalid_networks}

    analyze_networks(networks_to_analyze, network_polyploids)
