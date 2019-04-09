import unittest

import os

from eva_cttv_pipeline.evidence_string_generation import clinvar_to_evidence_strings
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT
from eva_cttv_pipeline.evidence_string_generation import trait
from tests.evidence_string_generation import test_clinvar
from tests.evidence_string_generation import config


def _get_mappings():
    efo_mapping_file = os.path.join(os.path.dirname(__file__), 'resources',
                                    'feb16_jul16_combined_trait_to_url.tsv')
    ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
    snp_2_gene_file = os.path.join(os.path.dirname(__file__), 'resources',
                                   config.snp_2_gene_file)

    mappings = clinvar_to_evidence_strings.get_mappings(efo_mapping_file, snp_2_gene_file)

    return mappings


MAPPINGS = _get_mappings()


class GetMappingsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mappings = MAPPINGS

    def test_efo_mapping(self):
        self.assertEqual(len(self.mappings.trait_2_efo), 5283)

        self.assertEqual(self.mappings.trait_2_efo["renal-hepatic-pancreatic dysplasia 2"][0],
                         ('http://www.orpha.net/ORDO/Orphanet_294415', None))
        self.assertEqual(self.mappings.trait_2_efo["frontotemporal dementia"][0],
                         ('http://purl.obolibrary.org/obo/HP_0000733', None))
        self.assertEqual(
            self.mappings.trait_2_efo["3 beta-hydroxysteroid dehydrogenase deficiency"][0],
            ('http://www.orpha.net/ORDO/Orphanet_90791', None))

        self.assertEqual(
            self.mappings.trait_2_efo["coronary artery disease/myocardial infarction"],
            [('http://www.ebi.ac.uk/efo/EFO_0000612', 'myocardial infarction'),
             ('http://www.ebi.ac.uk/efo/EFO_0001645', 'coronary heart disease')])

    def test_consequence_type_dict(self):
        self.assertEqual(len(self.mappings.consequence_type_dict), 34)

        self.assertTrue("rs121908485" in self.mappings.consequence_type_dict)
        self.assertTrue("rs121912888" in self.mappings.consequence_type_dict)
        self.assertTrue("rs137852558" in self.mappings.consequence_type_dict)
        self.assertTrue("rs137853008" in self.mappings.consequence_type_dict)

        self.assertFalse("rs0" in self.mappings.consequence_type_dict)
        self.assertFalse("rs5" in self.mappings.consequence_type_dict)
        self.assertFalse("rs9" in self.mappings.consequence_type_dict)


class CreateTraitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.trait = clinvar_to_evidence_strings.create_trait_list(
            ["Ciliary dyskinesia, primary, 7"], MAPPINGS.trait_2_efo, 9)[0]

    def test_clinvar_trait_list(self):
        self.assertEqual(self.trait.clinvar_name, 'ciliary dyskinesia, primary, 7')

    def test_efo_list(self):
        self.assertEqual(self.trait.ontology_id, 'http://www.ebi.ac.uk/efo/EFO_0003900')

    def test_multiple_mappings(self):
        trait1 = trait.Trait("barrett esophagus/esophageal adenocarcinoma",
                             "http://www.ebi.ac.uk/efo/EFO_0000478",
                             "esophageal adenocarcinoma", 1)
        trait2 = trait.Trait("barrett esophagus/esophageal adenocarcinoma",
                             "http://www.ebi.ac.uk/efo/EFO_0000280",
                             "Barrett's esophagus", 1)

        test_trait_list = clinvar_to_evidence_strings.create_trait_list(
            ["barrett esophagus/esophageal adenocarcinoma"], MAPPINGS.trait_2_efo, 1)

        self.assertEqual([trait1, trait2], test_trait_list)


    def test_return_none(self):
        none_trait = clinvar_to_evidence_strings.create_trait_list(["not a real trait"],
                                                                   MAPPINGS.trait_2_efo, 9)
        self.assertIsNone(none_trait)


class SkipRecordTest(unittest.TestCase):

    #clinvar_record, clinvar_record_measure, consequence_type, allele_origin,
                # allowed_clinical_significance, report

    def setUp(self):
        self.clinvar_record = test_clinvar.get_test_record()
        report = clinvar_to_evidence_strings.Report()
        consequence_type = CT.ConsequenceType("ENSG00000163646", CT.SoTerm("stop_gained"))
        # skip_record(clinvarRecord, cellbase_record, allowed_clinical_significance, counters)
        self.args = [self.clinvar_record, self.clinvar_record.measures[0], consequence_type,
                     "germline", ["not provided"], report]
        # allowed clin sig changed to just "non provided" to match that in the test record

    def test_return_true(self):
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_rs_is_none(self):
        self.clinvar_record.rs = None
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_con_type_is_none(self):
        self.clinvar_record.consequence_type = None
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))


class LoadEfoMappingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
        efo_file = \
            os.path.join(os.path.dirname(__file__), 'resources', 'feb16_jul16_combined_trait_to_url.tsv')

        cls.trait_2_efo = clinvar_to_evidence_strings.load_efo_mapping(efo_file)
        cls.trait_2_efo_w_ignore = clinvar_to_evidence_strings.load_efo_mapping(efo_file)

    def test_just_mapping_trait_2_efo(self):
        self.assertEqual(len(self.trait_2_efo), 5283)


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


class TestGetDefaultAllowedClincalSignificance(unittest.TestCase):
    def test_get_default_allowed_clincal_significance(self):
        correct_list = ['unknown', 'untested', 'non-pathogenic', 'probable-non-pathogenic',
         'probable-pathogenic', 'pathogenic', 'drug-response', 'drug response',
         'histocompatibility', 'other', 'benign', 'protective', 'not provided',
         'likely benign', 'confers sensitivity', 'uncertain significance',
         'likely pathogenic', 'conflicting data from submitters', 'risk factor',
         'association']
        self.assertEqual(clinvar_to_evidence_strings.get_default_allowed_clinical_significance(),
                         correct_list)


class TestConvertAlleleOrigins(unittest.TestCase):
    def test_just_germline(self):
        orig_allele_origins = ["germline"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual(["germline"], converted_allele_origins)

    def test_just_somatic(self):
        orig_allele_origins = ["somatic"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual(["somatic"], converted_allele_origins)

    def test_just_tested_inconclusive(self):
        orig_allele_origins = ["tested-inconclusive"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([], converted_allele_origins)

    def test_just_other_germline(self):
        orig_allele_origins_list = [["unknown"],
                                    ["inherited"],
                                    ["maternal"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
            self.assertListEqual(["germline"], converted_allele_origins)

    def test_nonsense(self):
        orig_allele_origins_list = [["fgdsgfgs"],
                                    ["notarealorigin"],
                                    ["134312432:dasdfd"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
            self.assertListEqual([], converted_allele_origins)
        orig_allele_origins = ["fgdsgfgs", "germline"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual(["germline"], converted_allele_origins)

    def test_mixed_germline(self):
        orig_allele_origins_list = [["germline", "de novo"],
                                    ["germline", "inherited", "not applicable"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
            self.assertListEqual(["germline"], converted_allele_origins)

    def test_duplicate(self):
        orig_allele_origins = ["germline", "germline"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
            orig_allele_origins)
        self.assertListEqual(["germline"], converted_allele_origins)
        orig_allele_origins = ["inherited", "inherited", "germline"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
            orig_allele_origins)
        self.assertListEqual(["germline"], converted_allele_origins)
        orig_allele_origins = ["somatic", "somatic", "somatic"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
            orig_allele_origins)
        self.assertListEqual(["somatic"], converted_allele_origins)


    def test_stringcase(self):
        orig_allele_origins_list = [["Germline"],
                               ["InHerIted"],
                               ["UNKNOWN"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
            self.assertListEqual(["germline"], converted_allele_origins)
        orig_allele_origins_list = [["Somatic"],
                                    ["SOMATIC"],
                                    ["sOMatIc"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
                orig_allele_origins)
            self.assertListEqual(["somatic"], converted_allele_origins)

    def test_mixed(self):
        orig_allele_origins_list = [["germline", "somatic"],
                                    ["somatic", "inherited", "not applicable"],
                                    ["somatic", "unknown"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
                orig_allele_origins)
            self.assertListEqual(["somatic", "germline"], converted_allele_origins)


class TestGetConsequenceTypes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_crm = test_clinvar.get_test_record().measures[0]
        cls.consequence_type_dict = CT.process_consequence_type_file(config.snp_2_gene_file)

    def test_get_consequence_types(self):
        test_consequence_type = CT.ConsequenceType("ENSG00000163646", CT.SoTerm("stop_gained"))

        self.assertEqual(clinvar_to_evidence_strings.get_consequence_types(
            self.test_crm,
            self.consequence_type_dict)[0],
            test_consequence_type)
        self.assertEqual(clinvar_to_evidence_strings.get_consequence_types(self.test_crm, {}),
                         [None])

