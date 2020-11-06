from sys import argv
import pathlib
import csv


def convert_value(sval, empty_default=0):
    try:
        return float(sval)
    except ValueError:
        return empty_default


def convert_name(original_name):
    name_parts = []
    for word in original_name.split('_'):
        if not word.endswith('.'):
            name_parts.append(word)
    return ' '.join(name_parts)


def extract_polyploids(poly_path):
    # PLOIDY_COLUMNS 0 - tax name, 2 - sure polyploid, 8 - sus poly
    definite_polyploids, suspect_polyploids = set(), set()
    with open(poly_path) as poly_csv:
        poly_reader = csv.reader(poly_csv)
        next(poly_reader)

        for row in poly_reader:
            name, definite, suspect = convert_name(row[0]), convert_value(row[2]), \
                                      convert_value(row[8])
            if definite:
                definite_polyploids.add(name)
            if suspect:
                suspect_polyploids.add(name)

    return definite_polyploids, suspect_polyploids


if __name__ == "__main__":
    polyploid_path = pathlib.Path(argv[1])
    definite_polyploids, suspect_polyploids = extract_polyploids(polyploid_path)
    print(len(definite_polyploids), len(suspect_polyploids))
