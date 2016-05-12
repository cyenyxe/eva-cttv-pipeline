from datetime import datetime
import json
import math
import os
import unittest

from eva_cttv_pipeline import clinvar_to_evidence_strings, config, consequence_type
from tests import test_clinvar_record


class SkipRecordTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.clinvar_record = test_clinvar_record.get_test_record()

    def setUp(self):
        report = clinvar_to_evidence_strings.Report()
        self.record = clinvar_to_evidence_strings.create_record({"reference": "A", "alternate": "T"}, None,
                                                                clin_sig="pathogenic",
                                                                clinvarRecord=self.clinvar_record,
                                                                con_type="transcript_ablation", rs="rs1")
        # skip_record(cellbase_record, record, allowed_clinical_significance, rcv_to_nsv, counters)
        self.args = [self.record, ["pathogenic", "likely pathogenic"],
                     {'RCV000138025': 'nsv869213', 'RCV000133922': 'nsv491994'}, report.counters]

    def test_return_false(self):
        self.assertFalse(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_not_in_allowed_clinical_significance(self):
        self.record.clin_sig = "unknown"
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_ref_eq_alt(self):
        self.record.cellbase_record = {"reference": "A", "alternate": "A"}
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_rs_is_none(self):
        self.record.rs = None
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_con_type_is_none(self):
        self.record.con_type = None
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

