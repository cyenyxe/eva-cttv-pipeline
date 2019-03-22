from functools import total_ordering, lru_cache
import json
import logging
import re
import requests

from eva_cttv_pipeline.trait_mapping.ols import get_ontology_label_from_ols, is_in_efo
from eva_cttv_pipeline.trait_mapping.ols import is_current_and_in_efo


logger = logging.getLogger(__package__)


class OntologyUri:
    db_to_uri_dict = {
        "orphanet": "http://www.orpha.net/ORDO/Orphanet_{}",
        "omim": "http://identifiers.org/omim/{}",
        "efo": "http://www.ebi.ac.uk/efo/EFO_{}",
        "mesh": "http://identifiers.org/mesh/{}",
        "medgen": "http://identifiers.org/medgen/{}",
        "hp": "http://purl.obolibrary.org/obo/HP_{}",
        "doid": "http://purl.obolibrary.org/obo/DOID_{}",
        "mondo": "http://purl.obolibrary.org/obo/MONDO_{}",
    }

    def __init__(self, id_, db):
        self.id_ = id_
        self.db = db
        self.uri = self.db_to_uri_dict[self.db.lower()].format(self.id_)

    def __str__(self):
        return self.uri


@total_ordering
class OxOMapping:
    """
    Individual mapping for an ontology ID mapped to one other ontology ID. An OxO result can consist
    of multiple mappings.
    """
    def __init__(self, label, curie, distance, query_id):
        self.label = label
        self.curie = curie
        self.db, self.id_ = curie.split(":")
        self.uri = OntologyUri(self.id_, self.db)
        self.distance = distance
        self.query_id = query_id
        self.in_efo = False
        # For non-EFO mappings, `is_current` property does not make sense and it not used
        self.is_current = False
        self.ontology_label = ""

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        return (self.label == other.label, self.db == other.db, self.id_ == other.id_,
                self.distance == other.distance, self.in_efo == other.in_efo,
                self.is_current == other.is_current, self.ontology_label == other.ontology_label)

    def __lt__(self, other):
        return ((other.distance, self.in_efo, self.is_current) <
                (self.distance, other.in_efo, other.is_current))

    def __str__(self):
        return "{}, {}, {}, {}".format(self.label, self.curie, self.distance, self.query_id)


class OxOResult:
    """
    A single result from querying OxO for one ID. A result can contain multiple mappings. A response
    from OxO can contain multiple results- one per queried ID.
    """
    def __init__(self, query_id, label, curie):
        self.query_id = query_id
        self.label = label
        self.curie = curie
        self.db, self.id_ = curie.split(":")
        self.uri = OntologyUri(self.id_, self.db)
        self.mapping_list = []

    def __str__(self):
        return "{}, {}, {}".format(self.query_id, self.label, self.curie)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        return (self.query_id == other.query_id, self.label == other.label,
                self.curie == other.curie, self.db == other.db, self.id_ == other.id_,
                self.uri == other.uri, self.mapping_list == other.mapping_list)


URI_DB_TO_DB_DICT = {
    "ordo": "Orphanet",
    "orphanet": "Orphanet",
    "omim": "OMIM",
    "efo": "EFO",
    "mesh": "MeSH",
    "hp": "HP",
    "doid": "DOID",
    "mondo": "MONDO",
}


NON_NUMERIC_RE = re.compile(r'[^\d]+')


@lru_cache(maxsize=16384)
def uri_to_oxo_format(uri: str) -> str:
    """
    Convert an ontology uri to a DB:ID format with which to query OxO

    :param uri: Ontology uri for a term
    :return: String in the format "DB:ID" with which to query OxO
    """
    if not any(x in uri.lower() for x in URI_DB_TO_DB_DICT.keys()):
        return None
    uri = uri.rstrip("/")
    uri_list = uri.split("/")
    if "identifiers.org" in uri:
        db = uri_list[-2]
        id_ = uri_list[-1]
    else:
        db, id_ = uri_list[-1].split("_")
    db = URI_DB_TO_DB_DICT[db.lower()]
    return "{}:{}".format(db, id_)


def uris_to_oxo_format(uri_set: set) -> list:
    """For each ontology uri in a set convert to the format of an ID suitable for querying OxO"""
    oxo_id_list = []
    for uri in uri_set:
        oxo_id = uri_to_oxo_format(uri)
        if oxo_id is not None:
            oxo_id_list.append(oxo_id)
    return oxo_id_list


def build_oxo_payload(id_list: list, target_list: list, distance: int) -> dict:
    """
    Build a dict containing the payload with which to make a POST request to OxO for finding xrefs
    for IDs in provided id_list, with the constraints provided in target_list and distance.

    :param id_list: List of IDs with which to find xrefs using OxO
    :param target_list: List of ontology datasources to include
    :param distance: Number of steps to take through xrefs to find mappings
    :return: dict containing payload to be used in POST request with OxO
    """
    payload = {}
    payload["ids"] = id_list
    payload["mappingTarget"] = target_list
    payload["distance"] = distance
    return payload


def oxo_query_helper(url: str, payload: dict) -> dict:
    """
    Make post request to OxO url using provided payload, returning json response, or None if there
    is an error in decoding.

    :param url: url to make request
    :param payload: Payload to use to make POST request
    :return: json response from OxO
    """
    try:
        json_response = requests.post(url, data=payload).json()
        return json_response
    except json.decoder.JSONDecodeError as e:
        return None


def oxo_request_retry_helper(retry_count: int, url: str, id_list: list, target_list: list,
                             distance: int) -> dict:
    """
    Make a number of attempts to query OxO for it to successfully return a non-None value,
    subsequently returning this value. Makes the number of tries specified in retry_count parameter.

    :param retry_count: Number of attempts to make
    :param url: String specifying the url to make a request.
    :param id_list: List of IDs with which to find xrefs using OxO
    :param target_list: List of ontology datasources to include
    :param distance: Number of steps to take through xrefs to find mappings
    :return: Returned value from OxO request.
    """
    payload = build_oxo_payload(id_list, target_list, distance)
    for retry_num in range(retry_count):
        return_value = oxo_query_helper(url, payload)
        if return_value is not None:
            return return_value
        logger.warning("attempt {}: failed running function oxo_query_helper with url {}".format(
            retry_num, url))
    logger.warning("error on last attempt, skipping")
    return None


def get_oxo_results_from_response(oxo_response: dict) -> list:
    """
    For a json(/dict) response from an OxO request, parse the data into a list of OxOResults

    :param oxo_response: Response from OxO request
    :return: List of OxOResults based upon the response from OxO
    """
    oxo_result_list = []
    results = oxo_response["_embedded"]["searchResults"]
    for result in results:
        if len(result["mappingResponseList"]) == 0:
            continue
        query_id = result["queryId"]
        label = result["label"]
        curie = result["curie"]
        oxo_result = OxOResult(query_id, label, curie)
        for mapping_response in result["mappingResponseList"]:
            mapping_label = mapping_response["label"]
            mapping_curie = mapping_response["curie"]
            mapping_distance = mapping_response["distance"]
            oxo_mapping = OxOMapping(mapping_label, mapping_curie, mapping_distance, query_id)

            uri = str(oxo_mapping.uri)

            ontology_label = get_ontology_label_from_ols(uri)
            if ontology_label is not None:
                oxo_mapping.ontology_label = ontology_label

            uri_is_current_and_in_efo = is_current_and_in_efo(uri)
            if not uri_is_current_and_in_efo:
                uri_is_in_efo = is_in_efo(uri)
                oxo_mapping.in_efo = uri_is_in_efo
            else:
                oxo_mapping.in_efo = uri_is_current_and_in_efo
                oxo_mapping.is_current = uri_is_current_and_in_efo

            oxo_result.mapping_list.append(oxo_mapping)

        oxo_result_list.append(oxo_result)

    return oxo_result_list


def get_oxo_results(id_list: list, target_list: list, distance: int) -> list:
    """
    Use list of ontology IDs, datasource targets and distance call function to query OxO and return
    a list of OxOResults.

    :param id_list: List of ontology IDs with which to find xrefs using OxO
    :param target_list: List of ontology datasources to include
    :param distance: Number of steps to take through xrefs to find mappings
    :return: List of OxOResults based upon results from request made to OxO
    """
    url = "https://www.ebi.ac.uk/spot/oxo/api/search?size=5000"
    oxo_response = oxo_request_retry_helper(4, url, id_list, target_list, distance)

    if oxo_response is None:
        return []

    if "_embedded" not in oxo_response:
        logger.warning("Cannot parse the response from OxO for the following identifiers:")
        logger.warning(','.join(id_list))
        return []

    return get_oxo_results_from_response(oxo_response)
