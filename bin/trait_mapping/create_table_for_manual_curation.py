#!/usr/bin/env python3

import argparse

from eva_cttv_pipeline.trait_mapping.ols import (
    get_ontology_label_from_ols, is_current_and_in_efo, is_in_efo,
)


def find_previous_mapping(trait_name, previous_mappings):
    if trait_name not in previous_mappings:
        return ''
    uri = previous_mappings[trait_name]
    label = get_ontology_label_from_ols(uri)
    uri_is_current_and_in_efo = is_current_and_in_efo(uri)
    uri_in_efo = is_in_efo(uri)
    if uri_in_efo:
        trait_status = 'EFO_CURRENT' if uri_is_current_and_in_efo else 'EFO_OBSOLETE'
    else:
        trait_status = 'NOT_CONTAINED'
    trait_string = '|'.join([uri, label, 'NOT_SPECIFIED', 'previously-used', trait_status])
    return trait_string


def find_exact_mapping(trait_name, mappings):
    for mapping in mappings:
        if mapping.lower().split('|')[1] == trait_name:
            return mapping
    return ''


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--traits-for-curation',
        help='Table with traits for which the pipeline failed to make a confident prediction')
    parser.add_argument(
        '--previous-mappings',
        help='Table with all mappings previously issued by EVA')
    parser.add_argument(
        '--final-table-for-curation',
        help='TSV to be loaded in Google Sheets for manual curation')
    args = parser.parse_args()
    outfile = open(args.final_table_for_curation, 'w')

    # Load all previous mappings
    previous_mappings = dict(l.rstrip().split('\t') for l in open(args.previous_mappings))

    # Process all mappings which require manual curation
    for line in open(args.traits_for_curation):
        fields = line.rstrip().split('\t')
        trait_name, trait_freq = fields[:2]
        mappings = fields[2:]
        previous_mapping = find_previous_mapping(trait_name, previous_mappings)
        exact_mapping = find_exact_mapping(trait_name, mappings)
        out_line = '\t'.join(
            [trait_name, trait_freq,
             # Mapping to use, if ready, comment, mapping URI, mapping label, whether exact, in EFO
             '', '', '', '', '', '', '',
             previous_mapping, exact_mapping] + mappings
        ) + '\n'
        outfile.write(out_line)
