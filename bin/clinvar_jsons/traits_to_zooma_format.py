import argparse
import gzip
import itertools
import json
from functools import lru_cache
import logging
import sys
from time import gmtime, strftime

import progressbar
import requests

from clinvar_jsons_shared_lib import clinvar_jsons, get_traits_from_json, has_allowed_clinical_significance


DATE = strftime("%d/%m/%y %H:%M", gmtime())


def main(args):
    with open(args.outfile_path, "wt") as outfile:
        outfile.write("STUDY\tBIOENTITY\tPROPERTY_TYPE\tPROPERTY_VALUE\tSEMANTIC_TAG\tANNOTATOR\tANNOTATION_DATE\n")
        line_count = file_len(args.infile_path)
        bar = progressbar.ProgressBar(max_value=line_count,
                                      widgets=[progressbar.AdaptiveETA(samples=1000)])
        is_zooma_mapping_dict = {}
        for clinvar_json in bar(clinvar_jsons(args.infile_path)):
            if has_allowed_clinical_significance(clinvar_json):
                is_zooma_mapping_dict = process_clinvar_json(clinvar_json, outfile, args.zooma_host,
                                                             args.filters, is_zooma_mapping_dict)


def check_if_zooma_mapping(trait_name, is_zooma_mapping_dict, zooma_host, filters):

    if trait_name in is_zooma_mapping_dict:
        is_zooma_mapping = is_zooma_mapping_dict[trait_name]
    else:
        zooma_uri_set = get_zooma_uris(trait_name, zooma_host, filters)
        is_zooma_mapping = (zooma_uri_set is not None
                            and len(zooma_uri_set) > 0)
        is_zooma_mapping_dict[trait_name] = is_zooma_mapping

    return is_zooma_mapping


def process_clinvar_json(clinvar_json, outfile, zooma_host, filters, is_zooma_mapping_dict):
    clinvar_acc = get_clinvar_accession(clinvar_json)
    variant_id_list = get_variant_ids(clinvar_json)

    trait_dict = {}
    trait_dict = get_traits_from_json(clinvar_json, trait_dict)

    for variant_id, trait_dict_item in itertools.product(variant_id_list,
                                                         trait_dict.items()):
        trait_name, trait_obj = trait_dict_item
        if trait_name == "not provided":
            continue

        is_zooma_mapping = check_if_zooma_mapping(trait_name, is_zooma_mapping_dict, zooma_host,
                                                  filters)
        if is_zooma_mapping:
            continue

        for xref in trait_obj.xref_set:
            if xref.status.lower() != "current" or xref.db.lower() not in OntologyUri.db_to_uri_dict:
                continue
            ontology_uri = OntologyUri(xref.id_, xref.db)

            write_zooma_record(clinvar_acc, variant_id, trait_name, ontology_uri, DATE, outfile)

    return is_zooma_mapping_dict


def write_zooma_record(clinvar_acc, variant_id, trait_name, ontology_uri, date, outfile):
    zooma_output_list = [clinvar_acc,
                         variant_id,
                         "disease",
                         trait_name,
                         str(ontology_uri),
                         "clinvar-xrefs",
                         date]
    outfile.write("\t".join(zooma_output_list) + "\n")


def get_variant_ids(clinvar_json):
    variant_id_list = []
    if "measureSet" in clinvar_json["clinvarSet"]['referenceClinVarAssertion']:
        measure_list = clinvar_json["clinvarSet"]['referenceClinVarAssertion']["measureSet"]["measure"]
    elif "measureSet" in clinvar_json["clinvarSet"]['referenceClinVarAssertion']["genotypeSet"]:
        measure_list = []
        for measure_set in clinvar_json["clinvarSet"]['referenceClinVarAssertion']["genotypeSet"]["measureSet"]:
            for measure in measure_set["measure"]:
                measure_list.append(measure)
    else:
        raise KeyError()

    for measure in measure_list:
        if "xref" in measure:
            for xref in measure["xref"]:
                if xref["db"].lower() == "dbsnp":
                    variant_id_list.append("rs{}".format(xref["id"]))
                elif xref["db"].lower() == "dbvar" and xref["id"].lower()[:3] in ("nsv", "esv"):
                    variant_id_list.append(xref["id"])
    return variant_id_list


def get_clinvar_accession(clinvar_json):
    return clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["clinVarAccession"]["acc"]


def get_zooma_uris(trait_name, zooma_host, filters):
    url = build_zooma_query(trait_name, filters, zooma_host)
    json_response = request_retry_helper(zooma_query_helper, 4, url)

    if json_response is None:
        return None

    uri_set = set()
    for result in json_response:
        if result["confidence"].lower() == "high":
            uri_set.update(result["semanticTags"])

    return uri_set


def build_zooma_query(trait_name, filters, zooma_host):
    url = "{}/spot/zooma/v2/api/services/annotate?propertyValue={}".format(zooma_host, trait_name)
    url_filters = [
                    "required:[{}]".format(filters["required"]),
                    "ontologies:[{}]".format(filters["ontologies"]),
                    "preferred:[{}]".format(filters["preferred"])
                  ]
    url += "&filter={}".format(",".join(url_filters))
    return url


def request_retry_helper(function, retry_count, url):
    for retry_num in range(retry_count):
        return_value = function(url)
        if return_value is not None:
            return return_value
        print("attempt {}: failed running function {} with url {}".format(retry_num, function, url))
    print("error on last attempt, skipping")
    return None


@lru_cache(maxsize=16384)
def zooma_query_helper(url):
    try:
        json_response = requests.get(url).json()
        return json_response
    except json.decoder.JSONDecodeError as e:
        return None


def open_file(file_path, mode):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def file_len(fname):
    with open_file(fname, "rt") as f:
        for i, l in enumerate(f):
            pass
    return i + 1


class OntologyUri:
    db_to_uri_dict = {
        "orphanet": "http://www.orpha.net/ORDO/Orphanet_{}",
        "omim": "http://identifiers.org/omim/{}",
        "efo": "http://www.ebi.ac.uk/efo/{}",
        "mesh": "http://identifiers.org/mesh/{}",
        "medgen": "http://identifiers.org/medgen/{}",
        "human phenotype ontology": "http://purl.obolibrary.org/obo/HP_{}",
        "mondo": "http://purl.obolibrary.org/obo/MONDO_{}",
    }

    def __init__(self, id_, db):
        self.id_ = id_
        self.db = db
        # This separate condition for "human phenotype ontology" is needed because these IDs are
        # prefixed with "HP:" which should be ignored when creating uri
        if self.db.lower() == "human phenotype ontology":
            logging.warning('Found <human phenotype ontology>, this should not really happen')
            self.uri = self.db_to_uri_dict[self.db.lower()].format(self.id_[3:])
        else:
            self.uri = self.db_to_uri_dict[self.db.lower()].format(self.id_)

    def __str__(self):
        return self.uri


def parse_args(argv):
    description = """
            Script for extracting the trait names of ClinVar records from a file with a list
            of CellBase, ClinVar JSONs, and the number of traits with this trait name.
            """
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i", dest="infile_path", required=True, help="Path to a file containing one CellBase ClinVar JSON per line.'")
    parser.add_argument("-o", dest="outfile_path", required=True, help="Path to file to output trait names in the zooma-accepted format")

    parser.add_argument("-n", dest="ontologies", default="efo,ordo,hp,mondo",
                        help="ontologies to use in query")
    parser.add_argument("-r", dest="required", default="cttv,eva-clinvar,gwas",
                        help="data sources to use in query.")
    parser.add_argument("-p", dest="preferred", default="eva-clinvar,cttv,gwas",
                        help="preference for data sources, with preferred data source first.")
    parser.add_argument("-z", dest="zooma_host", default="https://www.ebi.ac.uk",
                        help="the host to use for querying zooma")

    args = parser.parse_args(args=argv[1:])

    args.infile_path = args.infile_path
    args.outfile_path = args.outfile_path

    args.filters = {"ontologies": args.ontologies, "required": args.required,
                    "preferred": args.preferred}

    args.zooma_host = args.zooma_host

    return args


if __name__ == '__main__':
    main(parse_args(sys.argv))
