from sys import argv
import pathlib
import random
import os

import pandas as pd
import numpy as np
from Bio import Phylo

from polyploid_species import PolyploidDictionary
from duplicate_finder import compare_all_networks
from network_analysis import rank_graph, network_nodf, get_fractional_indices
from phylogenetic_tree import verify_in_tree


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


def permute_table(table, row_order, column_order):
    ordered = table.reindex(row_order, axis='index')
    ordered = ordered.reindex(column_order, axis='columns')
    return ordered


def analyze_networks(networks, polyploids, phylogenetic_tree, result_path, measure='pg'):
    for name, table in networks.items():
        net_poly = polyploids[name]

        r_plants, r_pols = rank_graph(table, measure_by=measure)
        polyploid_indices = get_fractional_indices(r_plants, net_poly)
        """try:
            plants_fixed = ['_'.join(k.split()).lower() for k in r_plants.keys()]
            missing_in_tree = verify_in_tree(plants_fixed, phylogenetic_tree)  # verify with ALLMB tree in browser
            print(f"\n{missing_in_tree}/{len(r_plants)} of network {name} not found in tree.\n")
        except ValueError:
            print(f"failed to get whole tree for {name}")"""

        ordered_table = permute_table(
            table,
            sorted(r_plants.keys(), key=r_plants.get, reverse=True),
            sorted(r_pols, key=r_pols.get, reverse=True)
        )
        nodf, nodf_contributions = network_nodf(ordered_table)

        with open(result_path.joinpath(name), mode='w') as fp:
            fp.write(f"{len(r_plants)} {len(r_pols)}\n")  # total number of plants and pollinators for ref
            fp.write(f"{len(net_poly) / len(r_plants)}\n\n")  # polyploid fraction of plants

            for poly_sp in net_poly:  # raw polyploid measure
                fp.write(f"{r_plants[poly_sp]}\n")
            fp.write(f"\n{np.mean([r_plants[k] for k in net_poly])}\n\n")  # mean of polyploid measures

            for poly_sp in net_poly:  # polyploid importance indices by measure (higher == better)
                fp.write(f"{polyploid_indices[poly_sp]}\n")
            fp.write(f"\n{np.mean([v for k, v in polyploid_indices.items()])}\n\n")  # mean importance

            plant_mean_nodf = np.mean([v for k, v in nodf_contributions.items()])
            poly_mean_nodf = np.mean([nodf_contributions[k] for k in net_poly])
            try:
                fp.write(f"{poly_mean_nodf / plant_mean_nodf}\n")  # polyploid nodf contribution
            except FloatingPointError:
                assert poly_mean_nodf == 0.0 and plant_mean_nodf == 0.0, (poly_mean_nodf, plant_mean_nodf)
                fp.write(f"{0.0}\n")  # in case of no nestedness at all
                # this only happens in ponisio_2017_20140101_1319, which has 3 plants and 3 pollinators


if __name__ == "__main__":
    np.seterr(all='raise')
    analysis_results_path = pathlib.Path(argv[1])
    polyploid_path = pathlib.Path(argv[2])
    tree_path = pathlib.Path(argv[3]).joinpath('ALLMB.tre')
    network_paths = [pathlib.Path(path_str) for path_str in argv[4:]]

    polydict = PolyploidDictionary(polyploid_path)
    phylo_tree = Phylo.read(tree_path, format='newick')

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

    print(f"analyzing {len(networks_to_analyze)} networks")
    analyze_networks(networks_to_analyze, network_polyploids, phylo_tree, analysis_results_path,
                     measure='pg')
