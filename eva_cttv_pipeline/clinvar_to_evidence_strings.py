import json
import sys
import urllib.error
import urllib.parse
import urllib.request
import codecs

import jsonschema
import xlrd

from eva_cttv_pipeline import efo_term, clinvar_record, consequence_type, config
from eva_cttv_pipeline import evidence_strings as ES


__author__ = 'Javier Lopez: javild@gmail.com'


def clinvar_to_evidence_strings(dir_out, allowed_clinical_significance=None, ignore_terms_file=None,
                                adapt_terms_file=None, efo_mapping_file=None, snp_2_gene_file=None,
                                variant_summary_file=None):

    allowed_clinical_significance = allowed_clinical_significance.split(',') if allowed_clinical_significance else \
        ['unknown', 'untested', 'non-pathogenic', 'probable-non-pathogenic',
         'probable-pathogenic', 'pathogenic', 'drug-response', 'drug response',
         'histocompatibility', 'other', 'benign', 'protective', 'not provided',
         'likely benign', 'confers sensitivity', 'uncertain significance',
         'likely pathogenic', 'conflicting data from submitters', 'risk factor',
         'association']

    trait_2_efo, unavailable_efo_dict = load_efo_mapping(efo_mapping_file, ignore_terms_file, adapt_terms_file)

    consequence_type_dict = consequence_type.process_consequence_type_file(snp_2_gene_file)
    rcv_to_rs, rcv_to_nsv = clinvar_record.get_rcv_to_rsnsv_mapping(variant_summary_file)

    nsv_list = []
    evidence_string_list = []
    n_processed_clinvar_records = 0
    n_pathogenic_no_rs = 0
    n_multiple_evidence_strings = 0
    n_multiple_allele_origin = 0
    n_germline_somatic = 0
    n_records_no_recognised_allele_origin = 0
    no_variant_to_ensg_mapping = 0
    n_more_than_one_efo_term = 0
    n_same_ref_alt = 0
    n_missed_strings_unmapped_traits = 0
    n_nsvs = 0
    n_valid_rs_and_nsv = 0
    n_nsv_skipped_clin_sig = 0
    n_nsv_skipped_wrong_ref_alt = 0
    unmapped_traits = {}
    n_unrecognised_allele_origin = {}
    ensembl_gene_id_uris = set()
    traits = set()
    unrecognised_clin_sigs = set()
    evidence_list = []  # To store Hellen Parkinson records of the form [RCV, rs, ClinVar trait, EFO url]
    record_counter = 0
    n_total_clinvar_records = 0
    skip = 0

    answer = urllib.request.urlopen('http://' + config.HOST + '/cellbase/webservices/rest/v3/hsapiens/feature/clinical/all?source=clinvar&skip=' + str(skip) + '&limit=' + str(config.BATCH_SIZE))
    reader = codecs.getreader("utf-8")
    curr_response = json.load(reader(answer))['response'][0]
    curr_result_list = curr_response['result']

    # A progress bar is initialized
    # widgets = ['Loading evidence strings: ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
    # pbar = progressbar.ProgressBar(widgets=widgets, maxval=curr_response['numTotalResults']).start()
    print(str(curr_response['numTotalResults']) + ' ClinVar records in total.')
    while len(curr_result_list) > 0:
        n_total_clinvar_records += len(curr_result_list)
        for record in curr_result_list:
            record_counter += 1
            n_ev_strings_per_record = 0
            clinvarRecord = clinvar_record.ClinvarRecord(record['clinvarSet'])
            clin_sig = clinvarRecord.clinical_significance.lower()
            n_nsvs += (clinvarRecord.get_nsv(rcv_to_nsv) is not None)

            if clin_sig not in allowed_clinical_significance:
                if clinvarRecord.get_nsv(rcv_to_nsv) is not None:
                    n_nsv_skipped_clin_sig += 1
                continue

            if record['reference'] == record['alternate']:
                n_same_ref_alt += 1
                if clinvarRecord.get_nsv(rcv_to_nsv) is not None:
                    n_nsv_skipped_wrong_ref_alt += 1
                continue

            nsv_list = append_nsv(nsv_list, clinvarRecord, rcv_to_nsv)
            rs = clinvarRecord.get_rs(rcv_to_rs)

            if rs is None:
                n_pathogenic_no_rs += 1
                continue

            consequenceType = clinvarRecord.get_main_consequence_types(consequence_type_dict, rcv_to_rs)
            # Mapping rs->Gene was found at Mick's file and therefore ensembl_gene_id will never be None

            if consequenceType is None:
                no_variant_to_ensg_mapping += 1
                continue

            for ensembl_gene_id in consequenceType.ensembl_gene_ids:

                rcv_to_gene_evidence_codes = ['http://identifiers.org/eco/cttv_mapping_pipeline']  # Evidence codes provided by Mick
                ensembl_gene_id_uri = 'http://identifiers.org/ensembl/' + ensembl_gene_id
                trait_refs_list = [['http://europepmc.org/abstract/MED/' + str(ref) for ref in refList] for refList in clinvarRecord.trait_pubmed_refs]
                observed_regs_list = ['http://europepmc.org/abstract/MED/' + str(ref) for ref in clinvarRecord.observed_pubmed_refs]
                measure_set_refs_list = ['http://europepmc.org/abstract/MED/' + str(ref) for ref in clinvarRecord.measure_set_pubmed_refs]
                for trait_counter, trait_list in enumerate(clinvarRecord.traits):
                    clinvar_trait_list, efo_list = map_efo(trait_2_efo, trait_list)
                    # Only ClinVar records associated to a trait with mapped EFO term will generate evidence_strings
                    if len(efo_list) > 0:
                        clinvar_record_allele_origins = clinvarRecord.allele_origins
                        n_multiple_allele_origin += (len(clinvar_record_allele_origins) > 1)
                        n_germline_somatic += (('germline' in clinvar_record_allele_origins) and (
                        'somatic' in clinvar_record_allele_origins))
                        n_records_no_recognised_allele_origin += (
                        ('germline' not in clinvar_record_allele_origins) and (
                        'somatic' not in clinvar_record_allele_origins))
                        for allele_origin_counter, alleleOrigin in enumerate(clinvar_record_allele_origins):
                            if alleleOrigin == 'germline':
                                evidence_string = ES.CTTVGeneticsEvidenceString(efo_list,
                                                                                clin_sig,
                                                                                clinvarRecord,
                                                                                consequenceType,
                                                                                ensembl_gene_id,
                                                                                ensembl_gene_id_uri,
                                                                                measure_set_refs_list,
                                                                                observed_regs_list,
                                                                                rcv_to_gene_evidence_codes,
                                                                                record,
                                                                                rs,
                                                                                trait_counter,
                                                                                trait_refs_list,
                                                                                unrecognised_clin_sigs)
                                n_ev_strings_per_record = add_evidence_string(clinvarRecord, evidence_string,
                                                                              evidence_string_list,
                                                                              n_ev_strings_per_record)
                                evidence_list.append(
                                    [clinvarRecord.accession, rs, ','.join(clinvar_trait_list),
                                     ','.join(efo_list)])
                                n_valid_rs_and_nsv += (clinvarRecord.get_nsv(rcv_to_nsv) is not None)
                                n_more_than_one_efo_term += (len(efo_list) > 1)
                                traits.update(set(efo_list))
                                ensembl_gene_id_uris.add(ensembl_gene_id_uri)
                            elif alleleOrigin == 'somatic':
                                evidence_string = ES.CTTVSomaticEvidenceString(efo_list,
                                                                               clin_sig,
                                                                               clinvarRecord,
                                                                               ensembl_gene_id,
                                                                               ensembl_gene_id_uri,
                                                                               measure_set_refs_list,
                                                                               observed_regs_list,
                                                                               trait_counter,
                                                                               trait_refs_list,
                                                                               unrecognised_clin_sigs,
                                                                               consequenceType)
                                n_ev_strings_per_record = add_evidence_string(clinvarRecord, evidence_string,
                                                                              evidence_string_list,
                                                                              n_ev_strings_per_record)
                                evidence_list.append(
                                    [clinvarRecord.accession, rs, ','.join(clinvar_trait_list),
                                     ','.join(efo_list)])
                                n_valid_rs_and_nsv += (clinvarRecord.get_nsv(rcv_to_nsv) is not None)
                                n_more_than_one_efo_term += (len(efo_list) > 1)
                                traits.update(set(efo_list))
                                ensembl_gene_id_uris.add(ensembl_gene_id_uri)
                            elif alleleOrigin not in n_unrecognised_allele_origin:
                                n_unrecognised_allele_origin[alleleOrigin] = 1
                            else:
                                n_unrecognised_allele_origin[alleleOrigin] += 1
                    else:
                        n_missed_strings_unmapped_traits += 1
                        if trait_list[0] in unmapped_traits:
                            unmapped_traits[trait_list[0]] += 1
                        else:
                            unmapped_traits[trait_list[0]] = 1

                if n_ev_strings_per_record > 0:
                    n_processed_clinvar_records += 1
                    if n_ev_strings_per_record > 1:
                        n_multiple_evidence_strings += 1

            # pbar.update(record_counter)
        skip += config.BATCH_SIZE

        answer = urllib.request.urlopen('http://' + config.HOST + '/cellbase/webservices/rest/v3/hsapiens/feature/clinical/all?source=clinvar&skip=' + str(skip) + '&limit=' + str(config.BATCH_SIZE))
        reader = codecs.getreader("utf-8")
        curr_response = json.load(reader(answer))['response'][0]
        curr_result_list = curr_response['result']
    # pbar.finish()

    write_output(dir_out, nsv_list, unmapped_traits, unavailable_efo_dict, evidence_string_list, evidence_list)

    output_report(n_total_clinvar_records, evidence_string_list, n_processed_clinvar_records, unrecognised_clin_sigs,
                  n_same_ref_alt, ensembl_gene_id_uris, traits, n_pathogenic_no_rs, n_multiple_evidence_strings,
                  n_germline_somatic, n_multiple_allele_origin, n_unrecognised_allele_origin,
                  no_variant_to_ensg_mapping, n_missed_strings_unmapped_traits, n_records_no_recognised_allele_origin,
                  n_more_than_one_efo_term, unavailable_efo_dict, n_valid_rs_and_nsv, n_nsvs, n_nsv_skipped_clin_sig,
                  n_nsv_skipped_wrong_ref_alt)


def write_output(dir_out, nsv_list, unmapped_traits, unavailable_efo_dict, evidence_string_list, evidence_list):
    write_string_list_to_file(nsv_list, dir_out + '/' + config.NSV_LIST_FILE)

    fdw = open(dir_out + '/' + config.UNMAPPED_TRAITS_FILE_NAME, 'w')  # Contains traits without a mapping in Gary's xls
    fdw.write('Trait\tCount\n')
    for trait_list in unmapped_traits:
        fdw.write(str(trait_list.encode('utf8')) + '\t' + str(unmapped_traits[trait_list]) + '\n')
    fdw.close()

    fdw = open(dir_out + '/' + config.UNAVAILABLE_EFO_FILE_NAME,
               'w')  # Contains urls provided by Gary which are not yet included within EFO
    fdw.write('Trait\tCount\n')
    for url in unavailable_efo_dict:
        fdw.write(url.encode('utf8') + '\t' + str(unavailable_efo_dict[url]) + '\n')
    fdw.close()

    fdw = open(dir_out + '/' + config.EVIDENCE_STRINGS_FILE_NAME, 'w')
    for evidence_string in evidence_string_list:
        fdw.write(json.dumps(evidence_string) + '\n')
    fdw.close()

    fdw = open(dir_out + '/' + config.EVIDENCE_RECORDS_FILE_NAME, 'w')
    for evidenceRecord in evidence_list:
        fdw.write('\t'.join(evidenceRecord) + '\n')
    fdw.close()


def output_report(n_total_clinvar_records, evidence_string_list, n_processed_clinvar_records, unrecognised_clin_sigs,
                  n_same_ref_alt, ensembl_gene_id_uris, traits, n_pathogenic_no_rs, n_multiple_evidence_strings,
                  n_germline_somatic, n_multiple_allele_origin, n_unrecognised_allele_origin,
                  no_variant_to_ensg_mapping, n_missed_strings_unmapped_traits, n_records_no_recognised_allele_origin,
                  n_more_than_one_efo_term, unavailable_efo_dict, n_valid_rs_and_nsv, n_nsvs, n_nsv_skipped_clin_sig,
                  n_nsv_skipped_wrong_ref_alt):

    report_strings = [
        str(n_total_clinvar_records) + ' ClinVar records in total',
        str(len(evidence_string_list)) + ' evidence string jsons generated',
        str(n_processed_clinvar_records) + ' ClinVar records generated at least one evidence string',
        str(len(unrecognised_clin_sigs)) + " Clinical significance string(s) not found among those described in ClinVar documentation:",
        str(unrecognised_clin_sigs),
        str(n_same_ref_alt) + ' ClinVar records with allowed clinical significance did present the same reference and alternate and were skipped',
        'Activities of those ClinVar records with unrecognized clinical significances were set to "unknown".',
        str(len(ensembl_gene_id_uris)) + ' distinct ensembl gene ids appear in generated evidence string json objects',
            str(len(traits)) + ' distinct trait names found to include in generated evidence string json objects',
        str(n_pathogenic_no_rs) + ' ClinVar records with allowed clinical significance DO NOT have an rs id',
        str(n_multiple_evidence_strings) + ' ClinVar records generated more than one evidence_string',
        str(n_germline_somatic) + ' ClinVar records with germline and somatic origins',
        str(n_multiple_allele_origin) + ' ClinVar records with more than one allele origin',
        'Number valid ClinVar records with unprocessed allele origins:'
    ]

    report_strings.extend([' ' + alleleOrigin + ': ' + str(n_unrecognised_allele_origin[alleleOrigin]) for alleleOrigin in n_unrecognised_allele_origin])

    report_strings.extend([
        str(no_variant_to_ensg_mapping) + ' ClinVar records with allowed clinical significance and valid rs id were skipped due to a lack of Variant->ENSG mapping.',
        str(n_missed_strings_unmapped_traits) + ' ClinVar records with allowed clinical significance, valid rs id and Variant->ENSG mapping were skipped due to a lack of EFO mapping (see ' + config.UNMAPPED_TRAITS_FILE_NAME + ').',
        str(n_records_no_recognised_allele_origin) + ' ClinVar records with allowed clinical significance, valid rs id, valid Variant->ENSG mapping and valid EFO mapping were skipped due to a lack of a valid alleleOrigin.',
        str(n_more_than_one_efo_term) + ' evidence strings with more than one trait mapped to EFO terms',
        str(len(unavailable_efo_dict)) + ' evidence strings were generated with traits without EFO correspondence',
        str(n_valid_rs_and_nsv) + ' evidence strings were generated from ClinVar records with rs and nsv ids',
        str(n_nsvs) + ' total nsvs found',
        str(n_nsv_skipped_clin_sig) + ' ClinVar nsvs were skipped because of a different clinical significance',
        str(n_nsv_skipped_wrong_ref_alt) + ' ClinVar nsvs were skipped because of same ref and alt'
    ])

    for line in report_strings:
        print(line)


def add_evidence_string(clinvarRecord, ev_string, evidence_string_list, n_evidence_strings_per_record):
    try:
        ev_string.validate()
        evidence_string_list.append(ev_string)
        n_evidence_strings_per_record += 1
    except jsonschema.exceptions.ValidationError as err:
        print('Error: evidence_string does not validate against schema.')
        print('ClinVar accession: ' + clinvarRecord.accession)
        print(err)
        print(json.dumps(ev_string))
        sys.exit(1)
    except efo_term.EFOTerm.IsObsoleteException as err:
        print('Error: obsolete EFO term.')
        print('Term: ' + ev_string.get_disease().efoid)
        print(err)
        print(json.dumps(ev_string))
        sys.exit(1)

    return n_evidence_strings_per_record


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
