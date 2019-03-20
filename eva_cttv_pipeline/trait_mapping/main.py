from collections import Counter
import csv
import logging
import progressbar

from eva_cttv_pipeline.trait_mapping.output import output_trait
from eva_cttv_pipeline.trait_mapping.oxo import get_oxo_results
from eva_cttv_pipeline.trait_mapping.oxo import uris_to_oxo_format
from eva_cttv_pipeline.trait_mapping.trait import Trait
from eva_cttv_pipeline.trait_mapping.trait_names_parsing import parse_trait_names
from eva_cttv_pipeline.trait_mapping.zooma import get_zooma_results


logger = logging.getLogger(__package__)


def get_uris_for_oxo(zooma_result_list: list) -> set:
    """
    For a list of Zooma mappings return a list of uris for the mappings in that list with a high
    confidence.

    :param zooma_result_list: List with elements of class ZoomaResult
    :return: set of uris from high confidence Zooma mappings, for which to query OxO
    """
    uri_set = set()
    for mapping in zooma_result_list:
        # Only use high confidence Zooma mappings for querying OxO
        if mapping.confidence.lower() == "high":
            uri_set.update([entry.uri for entry in mapping.mapping_list])
    return uri_set


def process_trait(trait: Trait, filters: dict, zooma_host: str, oxo_target_list: list,
                  oxo_distance: int) -> Trait:
    """
    Process a single trait. Find any mappings in Zooma. If there are no high confidence Zooma
    mappings that are in EFO then query OxO with any high confidence mappings not in EFO.

    :param trait: The trait to be processed.
    :param filters: A dictionary of filters to use for querying Zooma.
    :param zooma_host: A string with the hostname to use for querying Zooma
    :param oxo_target_list: A list of strings, each being an OxO ID for an ontology. Used to specify
                            which ontologies should be queried using OxO.
    :param oxo_distance: int specifying the maximum number of steps to use to query OxO. i.e. OxO's
                         "distance" parameter.
    :return: The original trait after querying Zooma and possibly OxO, with any results found.
    """
    trait.zooma_result_list = get_zooma_results(trait.name, filters, zooma_host)
    trait.process_zooma_results()
    if (trait.is_finished
            or len(trait.zooma_result_list) == 0
            or any([entry.is_current
                    for mapping in trait.zooma_result_list
                    for entry in mapping.mapping_list])):
        return trait
    uris_for_oxo_set = get_uris_for_oxo(trait.zooma_result_list)
    if len(uris_for_oxo_set) == 0:
        return trait
    oxo_input_id_list = uris_to_oxo_format(uris_for_oxo_set)
    trait.oxo_result_list = get_oxo_results(oxo_input_id_list, oxo_target_list, oxo_distance)

    if not trait.oxo_result_list:
        logger.warning('No OxO mapping for trait {}'.format(trait.name))

    trait.process_oxo_mappings()

    return trait


def main(input_filepath, output_mappings_filepath, output_curation_filepath, filters, zooma_host,
         oxo_target_list, oxo_distance, unattended):
    trait_names_list = parse_trait_names(input_filepath)
    trait_names_counter = Counter(trait_names_list)

    with open(output_mappings_filepath, "w", newline='') as mapping_file, \
            open(output_curation_filepath, "wt") as curation_file:
        mapping_writer = csv.writer(mapping_file, delimiter="\t")
        mapping_writer.writerow(["#clinvar_trait_name", "uri", "label"])
        curation_writer = csv.writer(curation_file, delimiter="\t")

        trait_names_iterator = trait_names_counter.items()
        if not unattended:
            trait_names_iterator = progressbar.ProgressBar(
                trait_names_iterator, max_value=len(trait_names_counter),
                widgets=[progressbar.AdaptiveETA(samples=1000)]
            )

        logger.info("Loaded {} trait names".format(len(trait_names_counter)))
        for i, (trait_name, freq) in enumerate(trait_names_iterator):
            trait = Trait(trait_name, freq)
            trait = process_trait(trait, filters, zooma_host, oxo_target_list,
                                  oxo_distance)
            output_trait(trait, mapping_writer, curation_writer)
            if unattended and i % 10000 == 0:
                logger.info("Processed {} records".format(i))
