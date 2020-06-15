import requests as req

"""
    Taxonomy resolution using the GBIF backbone database (https://www.gbif.org/), used to
    order mangal data as a true bipartite network in the excel file.

    This is done by getting the kingdom of a given species: if it's from Plantae, it's a 
    plant, if it's from Animalia, it's a pollinator.
    
    Note that GBIF requests only work for either a 'genus' query or a 'species' query of
    two words (genus species).
"""

GBIF_URL = "https://api.gbif.org/v1/species/match?rank={rank}&strict=false&name={species_name}"


def gbif_request(search_value, delimeter='_', only_genus=False):
    request = GBIF_URL.format(rank=("genus" if only_genus else "species"),
                              species_name=delimeter.join(search_value.split(' ')))
    response = req.get(request)

    if response.status_code != 200:
        raise ConnectionError("request '{0}' failed with error code {1}.".format(request, response.status_code))
    return response.json()


def get_kingdom(name):
    if len(name.split(' ')) < 2:
        taxonomy = gbif_request(name, only_genus=True)
    elif len(name.split(' ')) == 2:
        taxonomy = gbif_request(name)
    else:
        raise ValueError("Request taxonomy only with 'species' or 'genus species'")

    if 'kingdom' not in taxonomy.keys():
        raise KeyError("{} kingdom not found, with response: {}".format(name, taxonomy))
    return taxonomy['kingdom']

def test_get_kingdom():
    assert get_kingdom("Echium wildpretii") == get_kingdom("Echium") == "Plantae"


if __name__ == "__main__":
    test_get_kingdom()
