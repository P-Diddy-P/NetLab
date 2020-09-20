from sys import argv
import os
import pathlib

import pandas as pd


def convert_table(src_file, dest_file):
    raw_table = pd.read_excel(src_file)
    plants = raw_table.iloc[0:1, :]
    pollinators = raw_table.iloc[:, 0:2]

    pollinators_clean = []
    for row_id, row in pollinators.iterrows():
        pollinators_clean.append(str(row[0]) + ' ' + str(row[1]))
    pollinators_clean = pollinators_clean[2:]

    plants_clean = []
    for label, value in plants.iteritems():
        plants_clean.append(str(label) + ' ' + str(value[0]))
    plants_clean = plants_clean[3:]

    clean_table = raw_table.drop(index=[0, 1]).drop(
        columns=raw_table.columns[0:3]).fillna(0)
    clean_table.columns = plants_clean
    clean_table[''] = pollinators_clean
    final_table = clean_table.set_index('').T
    final_table.to_csv(dest_file)


def convert_directory(src_dir, dest_dir):
    for src_name in os.listdir(src_dir):
        convert_table(
            src_dir.joinpath(src_name),
            dest_dir.joinpath(src_name).with_suffix('.csv')
        )
        print(f"done converting {src_name}")


if __name__ == "__main__":
    convert_directory(argv[1], argv[2])
