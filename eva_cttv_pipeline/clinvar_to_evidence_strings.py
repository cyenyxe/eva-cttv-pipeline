import itertools
import json
import sys
from collections import defaultdict
from types import SimpleNamespace

import jsonschema
import xlrd

from eva_cttv_pipeline import cellbase_records, efo_term, clinvar_record, consequence_type, config, evidence_strings


__author__ = 'Javier Lopez: javild@gmail.com'


class Report:
    def __init__(self, unavailable_efo_dict=None):
        if unavailable_efo_dict is None:
            self.unavailable_efo_dict = {}
        else:
            self.unavailable_efo_dict = unavailable_efo_dict

        self.unrecognised_clin_sigs = set()
        self.ensembl_gene_id_uris = set()
        self.traits = set()
        self.n_unrecognised_allele_origin = defaultdict(int)
        self.nsv_list = []
        self.unmapped_traits = defaultdict(int)
        self.evidence_string_list = []
        self.evidence_list = []  # To store Helen Parkinson records of the form
        self.counters = self._get_counters()

    def __str__(self):

        report_strings = [
            str(self.counters["record_counter"]) + ' ClinVar records in total',
            str(len(self.evidence_string_list)) + ' evidence string jsons generated',
            str(self.counters["n_processed_clinvar_records"]) + ' ClinVar records generated at least one evidence string',
            str(len(self.unrecognised_clin_sigs)) + " Clinical significance string(s) not found among those described in ClinVar documentation:",
            str(self.unrecognised_clin_sigs),
            str(self.counters["n_same_ref_alt"]) + ' ClinVar records with allowed clinical significance did present the same reference and alternate and were skipped',
            'Activities of those ClinVar records with unrecognized clinical significances were set to "unknown".',
            str(len(self.ensembl_gene_id_uris)) + ' distinct ensembl gene ids appear in generated evidence string json objects',
                str(len(self.traits)) + ' distinct trait names found to include in generated evidence string json objects',
            str(self.counters["n_pathogenic_no_rs"]) + ' ClinVar records with allowed clinical significance DO NOT have an rs id',
            str(self.counters["n_multiple_evidence_strings"]) + ' ClinVar records generated more than one evidence_string',
            str(self.counters["n_germline_somatic"]) + ' ClinVar records with germline and somatic origins',
            str(self.counters["n_multiple_allele_origin"]) + ' ClinVar records with more than one allele origin',
            'Number valid ClinVar records with unprocessed allele origins:'
        ]

        report_strings.extend([' ' + alleleOrigin + ': ' + str(self.n_unrecognised_allele_origin[alleleOrigin]) for alleleOrigin in self.n_unrecognised_allele_origin])

        report_strings.extend([
            str(self.counters["no_variant_to_ensg_mapping"]) + ' ClinVar records with allowed clinical significance and valid rs id were skipped due to a lack of Variant->ENSG mapping.',
            str(self.counters["n_missed_strings_unmapped_traits"]) + ' ClinVar records with allowed clinical significance, valid rs id and Variant->ENSG mapping were skipped due to a lack of EFO mapping (see ' + config.UNMAPPED_TRAITS_FILE_NAME + ').',
            str(self.counters["n_records_no_recognised_allele_origin"]) + ' ClinVar records with allowed clinical significance, valid rs id, valid Variant->ENSG mapping and valid EFO mapping were skipped due to a lack of a valid alleleOrigin.',
            str(self.counters["n_more_than_one_efo_term"]) + ' evidence strings with more than one trait mapped to EFO terms',
            str(len(self.unavailable_efo_dict)) + ' evidence strings were generated with traits without EFO correspondence',
            str(self.counters["n_valid_rs_and_nsv"]) + ' evidence strings were generated from ClinVar records with rs and nsv ids',
            str(self.counters["n_nsvs"]) + ' total nsvs found',
            str(self.counters["n_nsv_skipped_clin_sig"]) + ' ClinVar nsvs were skipped because of a different clinical significance',
            str(self.counters["n_nsv_skipped_wrong_ref_alt"]) + ' ClinVar nsvs were skipped because of same ref and alt'
        ])

        return '\n'.join(report_strings)

    def add_evidence_string(self, ev_string):
        try:
            ev_string.validate()
            self.evidence_string_list.append(ev_string)
        except jsonschema.exceptions.ValidationError as err:
            print('Error: evidence_string does not validate against schema.')
            # print('ClinVar accession: ' + record.clinvarRecord.accession)
            print(err)
            print(json.dumps(ev_string))
            sys.exit(1)
        except efo_term.EFOTerm.IsObsoleteException as err:
            print('Error: obsolete EFO term.')
            print('Term: ' + ev_string.get_disease().efoid)
            print(err)
            print(json.dumps(ev_string))
            sys.exit(1)

    def write_output(self, dir_out):
        write_string_list_to_file(self.nsv_list, dir_out + '/' + config.NSV_LIST_FILE)

        fdw = open(dir_out + '/' + config.UNMAPPED_TRAITS_FILE_NAME, 'w')  # Contains traits without a mapping in Gary's xls
        fdw.write('Trait\tCount\n')
        for trait_list in self.unmapped_traits:
            fdw.write(str(trait_list.encode('utf8')) + '\t' + str(self.unmapped_traits[trait_list]) + '\n')
        fdw.close()

        fdw = open(dir_out + '/' + config.UNAVAILABLE_EFO_FILE_NAME,
                   'w')  # Contains urls provided by Gary which are not yet included within EFO
        fdw.write('Trait\tCount\n')
        for url in self.unavailable_efo_dict:
            fdw.write(url.encode('utf8') + '\t' + str(self.unavailable_efo_dict[url]) + '\n')
        fdw.close()

        fdw = open(dir_out + '/' + config.EVIDENCE_STRINGS_FILE_NAME, 'w')
        for evidence_string in self.evidence_string_list:
            fdw.write(json.dumps(evidence_string) + '\n')
        fdw.close()

        fdw = open(dir_out + '/' + config.EVIDENCE_RECORDS_FILE_NAME, 'w')
        for evidenceRecord in self.evidence_list:
            fdw.write('\t'.join(evidenceRecord) + '\n')
        fdw.close()

    @staticmethod
    def _get_counters():
        return {"n_processed_clinvar_records" : 0, "n_pathogenic_no_rs": 0, "n_multiple_evidence_strings": 0,
            "n_multiple_allele_origin": 0, "n_germline_somatic": 0, "n_records_no_recognised_allele_origin": 0,
            "no_variant_to_ensg_mapping": 0, "n_more_than_one_efo_term": 0, "n_same_ref_alt": 0,
            "n_missed_strings_unmapped_traits": 0, "n_nsvs": 0, "n_valid_rs_and_nsv": 0, "n_nsv_skipped_clin_sig": 0,
            "n_nsv_skipped_wrong_ref_alt": 0, "record_counter": 0, "n_total_clinvar_records": 0}


def launch_pipeline(dir_out, allowed_clinical_significance, ignore_terms_file, adapt_terms_file, efo_mapping_file,
                    snp_2_gene_file, variant_summary_file):

    allowed_clinical_significance = allowed_clinical_significance.split(',') if allowed_clinical_significance else \
        get_default_allowed_clincal_significance()

    mappings = get_mappings(efo_mapping_file, ignore_terms_file, adapt_terms_file, snp_2_gene_file,
                            variant_summary_file)

    report = clinvar_to_evidence_strings(allowed_clinical_significance, mappings)

    output(report, dir_out)


def output(report, dir_out):
    report.write_output(dir_out)
    print(report)


def clinvar_to_evidence_strings(allowed_clinical_significance, mappings):

    report = Report(mappings.unavailable_efo_dict)

    for cellbase_record in cellbase_records.get_records():

        record = create_record(cellbase_record, mappings)
        report.counters["record_counter"] += 1
        report.counters["n_nsvs"] += (record.clinvarRecord.get_nsv(mappings.rcv_to_nsv) is not None)
        append_nsv(report.nsv_list, record.clinvarRecord, mappings.rcv_to_nsv)

        if skip_record(record, allowed_clinical_significance, mappings.rcv_to_nsv, report.counters):
            continue

        report.counters["n_multiple_allele_origin"] += (len(record.clinvarRecord.allele_origins) > 1)
        report.counters["n_germline_somatic"] += (('germline' in record.clinvarRecord.allele_origins) and (
        'somatic' in record.clinvarRecord.allele_origins))
        report.counters["n_records_no_recognised_allele_origin"] += (
            ('germline' not in record.clinvarRecord.allele_origins) and
            ('somatic' not in record.clinvarRecord.allele_origins))

        traits = create_traits(record.clinvarRecord.traits, mappings.trait_2_efo, report)

        for ensembl_gene_id, trait, allele_origin in itertools.product(record.con_type.ensembl_gene_ids, traits,
                                                                       record.clinvarRecord.allele_origins):

            if allele_origin not in ('germline', 'somatic'):
                report.n_unrecognised_allele_origin[allele_origin] += 1
                continue
            else:
                if allele_origin == 'germline':
                    # todo look into if any of these arguments to ev strings can be removed and their use extracted out
                    evidence_string = evidence_strings.CTTVGeneticsEvidenceString(record,
                                                                                  report,
                                                                                  trait,
                                                                                  ensembl_gene_id)
                elif allele_origin == 'somatic':
                    evidence_string = evidence_strings.CTTVSomaticEvidenceString(record,
                                                                                 report,
                                                                                 trait,
                                                                                 ensembl_gene_id)
                report.add_evidence_string(record, evidence_string)
                report.evidence_list.append(
                    [record.clinvarRecord.accession, record.rs, ','.join(trait.clinvar_trait_list),
                     ','.join(trait.efo_list)])
                report.counters["n_valid_rs_and_nsv"] += (record.clinvarRecord.get_nsv(mappings.rcv_to_nsv) is not None)
                report.counters["n_more_than_one_efo_term"] += (len(trait.efo_list) > 1)
                report.traits.update(set(trait.efo_list))
                report.ensembl_gene_id_uris.add(evidence_strings.get_ensembl_gene_id_uri(ensembl_gene_id))

                record.n_ev_strings_per_record += 1

        if record.n_ev_strings_per_record > 0:
            report.counters["n_processed_clinvar_records"] += 1
            if record.n_ev_strings_per_record > 1:
                report.counters["n_multiple_evidence_strings"] += 1

    return report


def get_mappings(efo_mapping_file, ignore_terms_file, adapt_terms_file, snp_2_gene_file, variant_summary_file):
    mappings = SimpleNamespace()
    mappings.trait_2_efo, mappings.unavailable_efo_dict = load_efo_mapping(efo_mapping_file, ignore_terms_file,
                                                                           adapt_terms_file)

    mappings.consequence_type_dict = consequence_type.process_consequence_type_file(snp_2_gene_file)
    mappings.rcv_to_rs, mappings.rcv_to_nsv = clinvar_record.get_rcv_to_rsnsv_mapping(variant_summary_file)

    return mappings


def create_record(cellbase_record, mappings, **kwargs):
    record = SimpleNamespace()
    record.cellbase_record = cellbase_record
    record.n_ev_strings_per_record = 0

    record.clinvarRecord = clinvar_record.ClinvarRecord(cellbase_record['clinvarSet']) \
        if "clinvarRecord" not in kwargs else kwargs["clinvarRecord"]

    record.clin_sig = record.clinvarRecord.clinical_significance.lower() \
        if "clin_sig" not in kwargs else kwargs["clin_sig"]

    record.rs = record.clinvarRecord.get_rs(mappings.rcv_to_rs) \
        if "rs" not in kwargs else kwargs["rs"]

    record.con_type = record.clinvarRecord.get_main_consequence_types(mappings.consequence_type_dict, mappings.rcv_to_rs) \
        if "con_type" not in kwargs else kwargs["con_type"]
    # Mapping rs->Gene was found at Mick's file and therefore ensembl_gene_id will never be None

    record.trait_refs_list = [['http://europepmc.org/abstract/MED/' + str(ref) for ref in refList]
                              for refList in record.clinvarRecord.trait_pubmed_refs]
    record.observed_refs_list = ['http://europepmc.org/abstract/MED/' + str(ref)
                                 for ref in record.clinvarRecord.observed_pubmed_refs]
    record.measure_set_refs_list = ['http://europepmc.org/abstract/MED/' + str(ref)
                                    for ref in record.clinvarRecord.measure_set_pubmed_refs]

    return record


def skip_record(record, allowed_clinical_significance, rcv_to_nsv, counters):

    if record.clin_sig not in allowed_clinical_significance:
        if record.clinvarRecord.get_nsv(rcv_to_nsv) is not None:
            counters["n_nsv_skipped_clin_sig"] += 1
        return True

    if record.cellbase_record['reference'] == record.cellbase_record['alternate']:
        counters["n_same_ref_alt"] += 1
        if record.clinvarRecord.get_nsv(rcv_to_nsv) is not None:
            counters["n_nsv_skipped_wrong_ref_alt"] += 1
        return True

    if record.rs is None:
        counters["n_pathogenic_no_rs"] += 1
        return True

    if record.con_type is None:
        counters["no_variant_to_ensg_mapping"] += 1
        return True

    return False


def create_traits(clinvar_traits, trait_2_efo, report):
    traits = []
    for trait_counter, trait_list in enumerate(clinvar_traits):
        new_trait = create_trait(trait_counter, trait_list, trait_2_efo)
        if new_trait:
            traits.append(new_trait)
        else:
            report.counters["n_missed_strings_unmapped_traits"] += 1
            report.unmapped_traits[trait_list[0]] += 1
    return traits


def create_trait(trait_counter, trait_list, trait_2_efo):
    trait = SimpleNamespace()
    trait.trait_counter = trait_counter
    trait.clinvar_trait_list, trait.efo_list = map_efo(trait_2_efo, trait_list)
    # Only ClinVar records associated to a trait with mapped EFO term will generate evidence_strings
    if len(trait.efo_list) == 0:
        return None
    return trait


def write_string_list_to_file(string_list, filename):
    fd = open(filename, 'w')
    fd.write('\n'.join(string_list))
    fd.close()


def append_nsv(nsv_list, clinvarRecord, rcv_to_nsv):
    nsv = clinvarRecord.get_nsv(rcv_to_nsv)
    if nsv is not None:
        nsv_list.append(nsv)
    return nsv_list


def map_efo(trait_2_efo, trait_list):
    efo_list = []
    trait_list_to_return = []
    trait_string = trait_list[0].lower()
    if trait_string in trait_2_efo:
        for efo_trait in trait_2_efo[trait_string]:
            if efo_trait not in efo_list:  # First element in trait_list mus always be the "Preferred" trait name
                trait_list_to_return.append(trait_list[0])
                efo_list.append(efo_trait)
    else:
        for trait in trait_list[1:]:
            trait_string = trait.lower()
            if trait_string in trait_2_efo:
                for efo_trait in trait_2_efo[trait_string]:
                    if efo_trait not in efo_list:  # First element in trait_list mus always be the "Preferred" trait name
                        trait_list_to_return.append(trait)
                        efo_list.append(efo_trait)

    return trait_list_to_return, efo_list


def load_efo_mapping(efo_mapping_file, ignore_terms_file=None, adapt_terms_file=None):
    ignore_terms = get_terms_from_file(ignore_terms_file)
    adapt_terms = get_terms_from_file(adapt_terms_file)

    print('Loading phenotypes to EFO mapping...')
    efo_mapping_read_book = xlrd.open_workbook(efo_mapping_file, formatting_info=True)
    efo_mapping_read_sheet = efo_mapping_read_book.sheet_by_index(0)
    trait_2_efo = {}
    unavailable_efo = {}
    n_efo_mappings = 0
    for i in range(1, efo_mapping_read_sheet.nrows):
        if efo_mapping_read_sheet.cell_value(rowx=i, colx=1) != '':
            valid_efo, urls_to_adapt = get_urls(efo_mapping_read_sheet.cell_value(rowx=i, colx=1).split(', '), ignore_terms, adapt_terms)
            clinvar_trait = efo_mapping_read_sheet.cell_value(rowx=i, colx=0).lower()
            if len(valid_efo) > 0:
                trait_2_efo[clinvar_trait] = valid_efo
                n_efo_mappings += 1
            elif len(urls_to_adapt) > 0:
                trait_2_efo[clinvar_trait] = []
                for url in urls_to_adapt:
                    if url not in unavailable_efo:
                        unavailable_efo[url] = 1
                    else:
                        unavailable_efo[url] += 1
                    trait_2_efo[clinvar_trait].append(get_unmapped_url(url))

    print(str(n_efo_mappings) + ' EFO mappings loaded')
    print(str(len(unavailable_efo)) + ' urls without an actual valid EFO mapping')

    return trait_2_efo, unavailable_efo


class UnhandledUrlTypeException(Exception):
    pass


def get_unmapped_url(url):
    parts = url.split('/')
    if parts[-1].startswith("Orphanet_"):
        new_url = "http://purl.bioontology.org/ORDO/" + parts[-1]
    elif parts[-1].startswith("HP_"):
        new_url = "http://purl.bioontology.org/obo/" + parts[-1]
    else:
        raise UnhandledUrlTypeException("Error. Unhandled url type: " + url)

    return new_url


def get_urls(url_list, ignore_terms, adapt_terms):
    valid_efo = []
    urls_to_adapt = []
    for term in url_list:
        if term not in ignore_terms:
            if term in adapt_terms:
                urls_to_adapt.append(term)
            else:
                valid_efo.append(term)

    return valid_efo, urls_to_adapt


def get_terms_from_file(terms_file):
    if terms_file is not None:
        print('Loading list of terms...')
        fd = open(terms_file, 'r')
        terms_list = [line.rstrip() for line in fd]
        fd.close()
        print(str(len(terms_file)) + ' terms found at ' + terms_file)
    else:
        terms_list = []

    return terms_list


def get_default_allowed_clincal_significance():
    return ['unknown', 'untested', 'non-pathogenic', 'probable-non-pathogenic',
         'probable-pathogenic', 'pathogenic', 'drug-response', 'drug response',
         'histocompatibility', 'other', 'benign', 'protective', 'not provided',
         'likely benign', 'confers sensitivity', 'uncertain significance',
         'likely pathogenic', 'conflicting data from submitters', 'risk factor',
         'association']
