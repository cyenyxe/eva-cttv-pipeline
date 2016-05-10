from datetime import datetime
import json
import math
import os
import unittest

from eva_cttv_pipeline import clinvar_to_evidence_strings, config
from tests import test_clinvar_record


class CellbaseRecordsTest(unittest.TestCase):
    def test_get_curr_result_list(self):
        curr_result_list = clinvar_to_evidence_strings.get_curr_result_list(0)
        self.assertEqual(len(curr_result_list), config.BATCH_SIZE)

        curr_response = clinvar_to_evidence_strings.get_curr_response(0)
        curr_result_list = clinvar_to_evidence_strings.get_curr_result_list(curr_response['numTotalResults'])
        self.assertEqual(len(curr_result_list), 0)

        curr_result_list = clinvar_to_evidence_strings.get_curr_result_list(999999)
        self.assertEqual(len(curr_result_list), 0)

        curr_response = clinvar_to_evidence_strings.get_curr_response(0)
        len_to_expect = 20
        curr_result_list = clinvar_to_evidence_strings.get_curr_result_list(curr_response['numTotalResults'] - len_to_expect)
        self.assertEqual(len(curr_result_list), len_to_expect)

    # def test_get_curr_result_lists(self):
    #     curr_response = clinvar_to_evidence_strings.get_curr_response(0)
    #     num_total_results = curr_response['numTotalResults']
    #     list_counter = 0
    #     for list in clinvar_to_evidence_strings.get_curr_result_lists():
    #         list_counter += 1
    #     pred_num_lists = math.ceil(num_total_results / config.BATCH_SIZE)
    #     self.assertEqual(pred_num_lists, list_counter)


class SkipRecordTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.clinvar_record = test_clinvar_record.get_test_record()

    def setUp(self):
        counters = clinvar_to_evidence_strings.get_counters()
        self.args = [{"reference": "A", "alternate": "T"}, "pathogenic", ["pathogenic", "likely pathogenic"],
                     self.clinvar_record, {'RCV000138025': 'nsv869213', 'RCV000133922': 'nsv491994'}, "rs1",
                     "transcript_ablation", counters]

    def test_return_false(self):
        self.assertFalse(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_not_in_allowed_clinical_significance(self):
        self.args[1] = "unknown"
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_ref_eq_alt(self):
        self.args[0] = {"reference": "A", "alternate": "A"}
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_rs_is_none(self):
        self.args[5] = None
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_con_type_is_none(self):
        self.args[6] = None
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))



class LoadEfoMappingTest(unittest.TestCase):
    # TODO Make smaller files for testing, extracts from larger file. Ensure to create a smaller ignore file too, that matches a subset of the efo mapping file.
    @classmethod
    def setUpClass(cls):
        ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
        efo_file = os.path.join(os.path.dirname(__file__), 'resources', 'ClinVar_Traits_EFO_090915.xls')

        cls.trait_2_efo, cls.unavailable_efo = clinvar_to_evidence_strings.load_efo_mapping(efo_file)
        cls.trait_2_efo_w_ignore, cls.unavailable_efo_w_ignore = clinvar_to_evidence_strings.load_efo_mapping(efo_file, ignore_terms_file=ignore_file)

    def test_just_mapping_trait_2_efo(self):
        self.assertEqual(len(self.trait_2_efo), 3819)

    def test_w_ignore_trait_2_efo(self):
        self.assertEqual(len(self.trait_2_efo_w_ignore), 3528)


class GetUnmappedUrlTest(unittest.TestCase):
    def test_orphanet(self):
        url = "http://www.orpha.net/ORDO/Orphanet_2670"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url), "http://purl.bioontology.org/ORDO/Orphanet_2670")

    def test_hp(self):
        url = "http://purl.obolibrary.org/obo/HP_0000545"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url), "http://purl.bioontology.org/obo/HP_0000545")

    def test_bad_url(self):
        #TODO currently the function exits execution with bad urls
        pass


class GetTermsFromFileTest(unittest.TestCase):
    #TODO do the same for adapt terms file?
    @classmethod
    def setUpClass(cls):
        ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
        cls.ignore_terms = clinvar_to_evidence_strings.get_terms_from_file(ignore_file)

    def test_with_file(self):
        self.assertEqual(len(self.ignore_terms), 218)
        self.assertEqual(self.ignore_terms[0], "http://purl.obolibrary.org/obo/HP_0011677")
        self.assertEqual(self.ignore_terms[-1], "http://www.orpha.net/ORDO/Orphanet_120795")

    def test_no_file(self):
        self.assertEqual(clinvar_to_evidence_strings.get_terms_from_file(None), [])

