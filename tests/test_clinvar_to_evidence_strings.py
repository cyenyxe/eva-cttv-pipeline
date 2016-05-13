from datetime import datetime
import json
import math
import os
import unittest

from eva_cttv_pipeline import clinvar_to_evidence_strings, config, consequence_type, clinvar_record
from tests import test_clinvar_record, test_evidence_strings


def _get_mappings():
    efo_mapping_file = os.path.join(os.path.dirname(__file__), 'resources', 'ClinVar_Traits_EFO_090915.xls')
    ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
    snp_2_gene_file = os.path.join(os.path.dirname(__file__), 'resources', 'cttv012_snp2gene_20160222.tsv')
    variant_summary_file = os.path.join(os.path.dirname(__file__), 'resources', 'variant_summary_2015-05.txt')

    mappings = clinvar_to_evidence_strings.get_mappings(efo_mapping_file, ignore_file, None, snp_2_gene_file,
                                                        variant_summary_file)

    return mappings


MAPPINGS = _get_mappings()


class GetMappingsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mappings = MAPPINGS

    def test_efo_mapping(self):
        self.assertEqual(len(self.mappings.trait_2_efo), 3528)
        self.assertEqual(len(self.mappings.unavailable_efo_dict), 0)

        self.assertEqual(self.mappings.trait_2_efo["deafness, autosomal recessive 22"],
                         ['http://www.ebi.ac.uk/efo/EFO_0001063'])
        self.assertEqual(self.mappings.trait_2_efo["oculocutaneous albinism type 1b"],
                         ['http://www.orpha.net/ORDO/Orphanet_79434'])
        self.assertEqual(self.mappings.trait_2_efo["merosin deficient congenital muscular dystrophy"],
                         ['http://www.orpha.net/ORDO/Orphanet_258'])

    def test_consequence_type_dict(self):
        self.assertEqual(len(self.mappings.consequence_type_dict), 109367)

        self.assertTrue("rs724159824" in self.mappings.consequence_type_dict)
        self.assertTrue("rs34296458" in self.mappings.consequence_type_dict)
        self.assertTrue("rs199476100" in self.mappings.consequence_type_dict)
        self.assertTrue("rs80360485" in self.mappings.consequence_type_dict)

        self.assertFalse("rs0" in self.mappings.consequence_type_dict)
        self.assertFalse("rs5" in self.mappings.consequence_type_dict)
        self.assertFalse("rs9" in self.mappings.consequence_type_dict)

    def test_rcv_to_rs_nsv(self):
        self.assertEqual(len(self.mappings.rcv_to_rs), 134663)
        self.assertEqual(len(self.mappings.rcv_to_nsv), 14290)

        self.assertEqual(self.mappings.rcv_to_nsv["RCV000020147"], "nsv1067916")
        self.assertEqual(self.mappings.rcv_to_nsv["RCV000004182"], "nsv1067860")
        self.assertEqual(self.mappings.rcv_to_nsv["RCV000004183"], "nsv1067861")

        self.assertEqual(self.mappings.rcv_to_rs["RCV000061038"], "rs140870493")
        self.assertEqual(self.mappings.rcv_to_rs["RCV000038449"], "rs397517136")
        self.assertEqual(self.mappings.rcv_to_rs["RCV000126020"], "rs75686037")


class CreateRecordTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mappings = MAPPINGS
        cls.cellbase_record = test_evidence_strings._get_test_cellbase_record_gene()
        cls.clinvarRecord = clinvar_record.ClinvarRecord(cls.cellbase_record['clinvarSet'])
        cls.record = clinvar_to_evidence_strings.create_record(cls.cellbase_record, cls.mappings)

    def test_clinvar_record(self):
        self.assertEqual(self.record.clinvarRecord, self.cellbase_record['clinvarSet'])

    def test_clin_sig(self):
        self.assertEqual(self.record.clin_sig, self.clinvarRecord.clinical_significance.lower())
        self.assertEqual(self.record.clin_sig, "likely pathogenic")

    def test_rs(self):
        self.assertEqual(self.record.rs, self.clinvarRecord.get_rs(self.mappings.rcv_to_rs))
        self.assertEqual(self.record.rs, "rs515726230")

    def test_con_type(self):
        self.assertEqual(self.record.con_type, self.clinvarRecord.get_main_consequence_types(
            self.mappings.consequence_type_dict, self.mappings.rcv_to_rs))

    def test_trait_refs_list(self):
        trait_refs_list_t = [['http://europepmc.org/abstract/MED/' + str(ref) for ref in refList] for refList in
                             self.clinvarRecord.trait_pubmed_refs]
        self.assertEqual(self.record.trait_refs_list, trait_refs_list_t)

    def test_observed_refs_list(self):
        observed_refs_list_t = ['http://europepmc.org/abstract/MED/' + str(ref)
                                for ref in self.clinvarRecord.observed_pubmed_refs]
        self.assertEqual(self.record.observed_refs_list, observed_refs_list_t)

    def test_measure_set_refs_list(self):
        measure_set_refs_list_t = ['http://europepmc.org/abstract/MED/' + str(ref)
                                   for ref in self.clinvarRecord.observed_pubmed_refs]
        self.assertEqual(self.record.measure_set_refs_list, measure_set_refs_list_t)


class CreateTraitTest(unittest.TestCase):
    pass


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
