from sys import argv
import pathlib
import os

import pandas as pd


def iter_sources(sources):
    for src in sources:
        for file in os.listdir(src):
            if file.endswith('csv'):
                yield src.joinpath(file)


def extract_networks(source_directories):
    all_networks = dict()
    for net_file in iter_sources(source_directories):
        print(net_file.stem)
        print(pd.read_csv(net_file, index_col=0))


if __name__ == "__main__":
    extract_networks([pathlib.Path(s) for s in argv[1:]])
