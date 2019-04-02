#!/usr/bin/env python3

import sys
from eva_cttv_pipeline.evidence_string_generation import utilities, clinvar_to_evidence_strings


def main():
    parser = utilities.ArgParser(sys.argv)
    utilities.check_dir_exists_create(parser.out)
    clinvar_to_evidence_strings.launch_pipeline(
        parser.out, allowed_clinical_significance=parser.clinical_significance,
        efo_mapping_file=parser.efo_mapping_file, snp_2_gene_file=parser.snp_2_gene_file,
        json_file=parser.json_file, ot_schema=parser.ot_schema)


if __name__ == '__main__':
    main()
