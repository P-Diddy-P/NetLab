import pathlib
from sys import argv


def extract_networks(source_directories):
    for src in source_directories:
        print(src)


if __name__ == "__main__":
    extract_networks([pathlib.Path(s) for s in argv[1:]])
