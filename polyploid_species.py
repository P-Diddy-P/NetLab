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


class PolyploidDictionary:
    def __init__(self, poly_path):
        # PLOIDY_COLUMNS 0 - tax name, 2 - sure polyploid, 8 - sus poly
        self.definite_polyploids, self.suspect_polyploids = set(), set()
        with open(poly_path) as poly_csv:
            poly_reader = csv.reader(poly_csv)
            next(poly_reader)

            for row in poly_reader:
                name, definite, suspect = convert_name(row[0]), convert_value(row[2]), \
                                      convert_value(row[8])
                if definite:
                    self.definite_polyploids.add(name)
                if suspect:
                    self.suspect_polyploids.add(name)

    def test_ploidy(self, tax_name, annotate=False):
        ploidy = 0
        if tax_name in self.definite_polyploids:
            ploidy = 2
        if tax_name in self.suspect_polyploids:
            ploidy = 1

        if annotate:
            ploidy = ploidy * '*'
        return ploidy

    @property
    def definite(self):
        return self.definite_polyploids

    @property
    def suspect(self):
        return self.suspect_polyploids


if __name__ == "__main__":
    polyploid_path = pathlib.Path(argv[1])
    polydict = PolyploidDictionary(polyploid_path)
    print(len(polydict.definite), len(polydict.suspect))
