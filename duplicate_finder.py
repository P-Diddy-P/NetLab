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


if __name__ == "__main__":
    # TODO currently doesn't parse IWDB networks, as they have no set format.
    # TODO add IWDB (argv[1]) once it's directory is in order.
    networks = extract_networks([pathlib.Path(s) for s in argv[2:]])
    network_keys = list(networks.keys())

    for i in range(len(network_keys)):
        refnet = networks[network_keys[i]]
        refplants = refnet.index  # add preprocessing for vertex lists
        refpols = refnet.columns

        for j in range(i + 1, len(network_keys)):
            pass
