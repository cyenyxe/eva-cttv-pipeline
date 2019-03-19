import argparse
import json
import requests
import sys


class Trait:

    def __init__(self, name, xref_string, frequency):
        self.name = name
        self.xref_string = xref_string
        self.frequency = frequency

    def __str__(self):
        return "{}\t{}\t{}".format(self.name, self.xref_string, self.frequency)


class OntologyMapping:

    def __init__(self, uris, zooma_label, confidence, source):
        self.uris = uris
        self.zooma_label = zooma_label
        self.labels = [self.zooma_label for _ in range(len(uris))]
        self.confidence = confidence
        self.source = source
        self.ols_label = ["0" for _ in range(len(uris))]

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}".format("|".join(self.labels), "|".join(self.ols_label), "|".join(self.uris), self.confidence, self.source)


def main():
    parser = ArgParser(sys.argv)

    # Read trait names from input file
    traits = read_traits(parser.input_filepath)

    with open(parser.output_filepath, "wt") as output_file:
        output_file.write("#{}\n".format(parser.filters))
        for trait in traits:
            output_file.write(str(trait))
            mappings = get_ontology_mappings(trait, parser.filters, parser.zooma_host)
            if mappings is not None:
                for mapping in mappings:
                    output_file.write("\t" + str(mapping))
            else:
                output_file.write("\tZOOMA_MAPPING_FAILED")
            output_file.write("\n")


def read_traits(filepath):
    traits = []
    with open(filepath) as f:
        for line in f:
            line_list = line.rstrip().split("\t")
            name = line_list[0]
            xref_string = line_list[1]
            frequency = line_list[2]
            traits.append(Trait(name, xref_string, frequency))
    return traits


def request_retry_helper(function, retry_count, url):
    for retry_num in range(retry_count):
        return_value = function(url)
        if return_value is not None:
            return return_value
        print("attempt {}: failed running function {} with url {}".format(retry_num, function, url))
    print("error on last attempt, skipping")
    return None


def zooma_query_helper(url):
    try:
        json_response_1 = requests.get(url).json()
        return json_response_1
    except json.decoder.JSONDecodeError as e:
        return None


def ols_query_helper(url):
    try:
        json_response = requests.get(url).json()
        for term in json_response["_embedded"]["terms"]:
            if term["is_defining_ontology"]:
                return term["label"]
    except:
        return None


def get_ontology_mappings(trait, filters, zooma_host):
    '''
    First get the URI, label from a selected source, confidence and source:
    http://snarf.ebi.ac.uk:8580/spot/zooma/v2/api/services/annotate?propertyValue=intellectual+disability
    Then the ontology label to replace the label from a source:
    https://www.ebi.ac.uk/ols/api/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0003847
    '''
    url = build_zooma_query(trait.name, filters, zooma_host)
    json_response_1 = request_retry_helper(zooma_query_helper, 4, url)

    if json_response_1 is None:
        return None

    mappings = get_mappings_for_trait(json_response_1)

    for mapping in mappings:
        for idx, uri in enumerate(mapping.uris):
            label = get_ontology_label_from_ols(uri)
            # If no label is returned (shouldn't really happen) keep the existing one
            if label is not None:
                mapping.labels[idx] = label
                mapping.ols_label[idx] = "1"
            else:
                print(
                    "Couldn't retrieve ontology label from OLS for trait '{}', will use the one from Zooma".format(
                        trait.name))

    return mappings


def build_zooma_query(trait_name, filters, zooma_host):
    url = "{}/spot/zooma/v2/api/services/annotate?propertyValue={}".format(zooma_host, trait_name)
    url_filters = [
                    "required:[{}]".format(filters["required"]),
                    "ontologies:[{}]".format(filters["ontologies"]),
                    "preferred:[{}]".format(filters["preferred"])
                  ]
    url += "&filter={}".format(",".join(url_filters))
    return url


def get_mappings_for_trait(zooma_response):
    mappings = []
    for result in zooma_response:
        uris = result["semanticTags"]
        zooma_label = result["annotatedProperty"]["propertyValue"]
        confidence = result["confidence"]
        source_name = result["derivedFrom"]["provenance"]["source"]["name"]
        mappings.append(OntologyMapping(uris, zooma_label, confidence, source_name))
    return mappings


def get_ontology_label_from_ols(uri_mapping):
    url = build_ols_query(uri_mapping)
    json_response_1 = request_retry_helper(ols_query_helper, 4, url)
    return json_response_1


def build_ols_query(ontology_uri):
    url = "https://www.ebi.ac.uk/ols/api/terms?iri={}".format(ontology_uri)
    return url


class ArgParser:

    def __init__(self, argv):
        description = """
                Script for running terms through Zooma, retrieving mapped uri, label from OLS,
                confidence of the mapping, and source of the mapping.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="input_filepath", required=True, help="path to input file, with trait names in first column, number of variants the trait name appears in in the second column. delimeted using tab")
        parser.add_argument("-o", dest="output_filepath", required=True, help="path to output file (not just the directory). outputs a file with a header (line starting with \"#\") which shows the filters used. then the first column is trait name, then number of variants for the trait, then zooma label, uri(s), confidence, source. these zooma columns repeat when there are multiple mappings.")
        parser.add_argument("-n", dest="ontologies", default="efo,ordo,hp,mondo", help="ontologies to use in query")
        parser.add_argument("-r", dest="required", default="cttv,eva-clinvar,gwas", help="data sources to use in query.")
        parser.add_argument("-p", dest="preferred", default="eva-clinvar,cttv,gwas", help="preference for data sources, with preferred data source first.")
        parser.add_argument("-z", dest="zooma_host", default="http://snarf.ebi.ac.uk:8580", help="the host to use for querying zooma")  # alternate to default is https://www.ebi.ac.uk

        args = parser.parse_args(args=argv[1:])

        self.input_filepath = args.input_filepath
        self.output_filepath = args.output_filepath

        self.filters = {"ontologies": args.ontologies, "required": args.required, "preferred": args.preferred}

        self.zooma_host = args.zooma_host


if __name__ == "__main__":
    main()
