#!/usr/bin/env python3

import argparse
import re
import requests

# Name of ontology in OLS url, e. g. https://www.ebi.ac.uk/ols/ontologies/ordo/terms?iri=...
ontology_to_ols = {
    'HP': 'hp',
    'MONDO': 'mondo',
    'Orphanet': 'ordo',
}

# OLS url to query for a term details
ols_url_template = 'https://www.ebi.ac.uk/ols/api/ontologies/{ontology}/terms?iri={term}'

# List of fields in the current version of Webulous submittion template
webulous_fields = [
    'disease', 'child_of', 'definition', 'synonyms', 'located_in_organ', 'located_in_cell',
    'biological_process', 'msh_def_cite', 'ncit_def_cite', 'snomedct_def_cite', 'icd9_def_cite',
    'icd10_def_cite', 'omim_def_cite', 'doid_def_cite', 'meddra_def_cite', 'umls_def_cite',
    'wikipedia_def_cite', 'comments', 'ordo_def_cite', 'definition_editor', 'definition_citation',
    'mondo_def_cite'
]
webulous_format_string = '\t'.join('{' + f + '}' for f in webulous_fields) + '\n'

# String to join multiple values in a Webulous template
webulous_joiner = ' || '


def get_parent_terms(url):
    return [term['label'] for term in requests.get(url).json()['_embedded']['terms']]


def get_ols_details(ontology, term):
    """Queries OLS and returns the details necessary for """
    url = ols_url_template.format(ontology=ontology, term=term)
    data = requests.get(url).json()['_embedded']['terms'][0]
    label = data['label']
    parents = get_parent_terms(data['_links']['parents']['href'])

    # Definition: first, try to get from annotation field
    definition = data['annotation'].get('definition', [''])[0]
    # If no luck, simply use a description
    if not definition and data['description']:
        definition = data['description'][0]

    synonyms = data['synonyms'] or []
    xrefs = {x.split(':')[0]: x for x in data['annotation'].get('database_cross_reference', [])}
    # If a term comes from either Orphanet or MONDO, we need to add these as xrefs as well
    # (since they won't be present in the normal list of xrefs).
    if ontology == 'mondo':
        xrefs['MONDO'] = term.split('/')[-1]
    elif ontology == 'ordo':
        xrefs['Orphanet'] = term.split('/')[-1].replace('_', ':')
    return label, parents, definition, synonyms, xrefs


def format_output_string(ontology, term):
    label, parents, definition, synonyms, xrefs = get_ols_details(ontology, term)

    return webulous_format_string.format(
        disease=label,
        child_of=webulous_joiner.join(parents),
        definition=definition,
        synonyms=webulous_joiner.join(synonyms),
        located_in_organ='',
        located_in_cell='',
        biological_process='',
        msh_def_cite=xrefs.get('MESH', ''),
        ncit_def_cite=xrefs.get('NCIT', ''),
        snomedct_def_cite='',
        icd9_def_cite='',
        icd10_def_cite='',
        omim_def_cite=xrefs.get('OMIM', ''),
        doid_def_cite=xrefs.get('DOID', ''),
        meddra_def_cite='',
        umls_def_cite=xrefs.get('UMLS', ''),
        wikipedia_def_cite='',
        comments='',
        ordo_def_cite=xrefs.get('Orphanet', ''),
        definition_editor='',
        definition_citation='',
        mondo_def_cite=xrefs.get('MONDO', ''),
    )


def create_efo_table(input_file_path, output_file_path):
    with open(input_file_path) as infile, open(output_file_path, 'w') as outfile:
        for line in infile:
            term = line.rstrip()
            print('Processing ' + term)
            ontology = ontology_to_ols[re.split('[:_]', term.split('/')[-1])[0]]
            result = format_output_string(ontology, term)
            outfile.write(result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-mappings', required=True,
                        help='Input file with ontology mappings')
    parser.add_argument('-o', '--output', required=True,
                        help='Output table for EFO import')
    args = parser.parse_args()
    create_efo_table(args.input_mappings, args.output)
