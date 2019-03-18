from enum import Enum
from functools import total_ordering, lru_cache
import json
import logging
import requests

from eva_cttv_pipeline.trait_mapping.ols import get_ontology_label_from_ols, \
    is_current_and_in_efo, is_in_efo
from eva_cttv_pipeline.trait_mapping.utils import request_retry_helper


@total_ordering
class ZoomaConfidence(Enum):
    """Enum to represent the confidence of a mapping in Zooma."""
    LOW = 1
    MEDIUM = 2
    GOOD = 3
    HIGH = 4

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        return self.value == other.value

    def __lt__(self, other):
        return self.value < other.value

    def __str__(self):
        return self.name


@total_ordering
class ZoomaMapping:
    """Representation of one ontology term in a mapping in Zooma."""
    def __init__(self, uri, confidence, source):
        self.uri = uri
        self.confidence = ZoomaConfidence[confidence.upper()]
        self.source = source
        self.ontology_label = ""
        self.in_efo = False
        self.is_current = False

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        if (self.uri != other.uri or self.confidence != other.confidence
                or self.ontology_label != other.ontology_label or self.in_efo != other.in_efo
                or self.is_current != other.is_current):
            return False
        return True

    def __lt__(self, other):
        return ((self.confidence, self.in_efo, self.is_current) <
                (other.confidence, other.in_efo, other.is_current))


class ZoomaResult:
    """
    A mapping in Zooma from one term, which can contain multiple ontology IDs mapped to. One
    term can be mapped to multiple mappings.
    """
    def __init__(self, uri_list, zooma_label, confidence, source):
        self.uri_list = uri_list
        self.zooma_label = zooma_label
        self.confidence = confidence
        self.source = source
        self.mapping_list = []
        for uri in uri_list:
            self.mapping_list.append(ZoomaMapping(uri, confidence, source))

    def __str__(self):
        return "{}, {}, {}, {}".format(self.zooma_label, self.confidence, self.source,
                                       self.mapping_list)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        return (self.uri_list == other.uri_list, self.zooma_label == other.zooma_label,
                self.confidence == other.confidence, self.source == other.source,
                self.mapping_list == other.mapping_list)


@lru_cache(maxsize=16384)
def zooma_query_helper(url: str) -> dict:
    """
    Make a get request to provided url and return the response, assumed to be a json response, in
    a dict.

    :param url: String of Zooma url used to make a request
    :return: Zooma response in a dict
    """
    try:
        json_response_1 = requests.get(url).json()
        return json_response_1
    except json.decoder.JSONDecodeError as e:
        return None


def get_zooma_results(trait_name: str, filters: dict, zooma_host: str) -> list:
    """
    Given a trait name, Zooma filters in a dict and a hostname to use, query Zooma and return a list
    of Zooma mappings for this trait.

    First get the URI, label from a selected source, confidence and source:
    http://snarf.ebi.ac.uk:8580/spot/zooma/v2/api/services/annotate?propertyValue=intellectual+disability
    Then the ontology label to replace the label from a source:
    https://www.ebi.ac.uk/ols/api/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0003847

    :param trait_name: A string containing a trait name from a ClinVar record.
    :param filters: A dictionary containing filters used when querying OxO
    :param zooma_host: Hostname of a Zooma instance to query.
    :return: List of ZoomaResults
    """
    url = build_zooma_query(trait_name, filters, zooma_host)
    zooma_response_list = request_retry_helper(zooma_query_helper, 4, url)

    if zooma_response_list is None:
        return []

    zooma_result_list = get_zooma_results_for_trait(zooma_response_list)

    for zooma_result in zooma_result_list:
        for zooma_mapping in zooma_result.mapping_list:
            label = get_ontology_label_from_ols(zooma_mapping.uri)
            # If no label is returned (shouldn't really happen) keep the existing one
            if label is not None:
                zooma_mapping.ontology_label = label
            else:
                logging.warning("Couldn't retrieve ontology label from OLS for trait '{}'".format(trait_name))

            uri_is_current_and_in_efo = is_current_and_in_efo(zooma_mapping.uri)
            if not uri_is_current_and_in_efo:
                uri_is_in_efo = is_in_efo(zooma_mapping.uri)
                zooma_mapping.in_efo = uri_is_in_efo
            else:
                zooma_mapping.in_efo = uri_is_current_and_in_efo
                zooma_mapping.is_current = uri_is_current_and_in_efo

    return zooma_result_list


def build_zooma_query(trait_name: str, filters: dict, zooma_host: str) -> str:
    """
    Given a trait name, filters and hostname, create a url with which to query Zooma. Return this
    url.

    :param trait_name: A string containing a trait name from a ClinVar record.
    :param filters: A dictionary containing filters used when querying OxO
    :param zooma_host: Hostname of a Zooma instance to query.
    :return: String of a url which can be requested
    """
    url = "{}/spot/zooma/v2/api/services/annotate?propertyValue={}".format(zooma_host, trait_name)
    url_filters = [
                    "required:[{}]".format(filters["required"]),
                    "ontologies:[{}]".format(filters["ontologies"]),
                    "preferred:[{}]".format(filters["preferred"])
                  ]
    url += "&filter={}".format(",".join(url_filters))
    return url


def get_zooma_results_for_trait(zooma_response_list: list) -> list:
    """
    Given a response from a Zooma request return ZoomaResults based upon the data in that request.

    :param zooma_response_list: A json (dict) response from a Zooma request.
    :return: List of ZoomaResulst in the Zooma response.
    """
    result_list = []
    for response in zooma_response_list:
        # uri_list = ",".join(result["semanticTags"])
        uris = response["semanticTags"]
        zooma_label = response["annotatedProperty"]["propertyValue"]
        confidence = response["confidence"]
        source_name = response["derivedFrom"]["provenance"]["source"]["name"]
        result_list.append(ZoomaResult(uris, zooma_label, confidence, source_name))
    return result_list
