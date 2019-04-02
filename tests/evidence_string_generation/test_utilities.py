import os
import unittest

from eva_cttv_pipeline.evidence_string_generation import utilities


class GetResourceFileTest(unittest.TestCase):
    def test_get_resource_file_existent(self):
        self.assertTrue(utilities.get_resource_file("eva_cttv_pipeline", "resources/json_schema"))

    def test_get_resource_file_nonexistent(self):
        self.assertEqual(utilities.get_resource_file("not_a_real_package_39146", "not_a_real_file"),
                         None)


class ArgParserTest(unittest.TestCase):
    clin_sig = 'pathogenic,likely pathogenic'
    ignore = '/path/to/ignore/file'
    out = '/path/to/out/file'
    efo_map_file = '/path/to/efo/file'
    snp_2_gene_file = '/path/to/snp/to/gene/file'
    ot_schema_path = '/path/to/ot/schema/json'

    @classmethod
    def setUpClass(cls):
        argv = ['clinvar_to_evidence_strings.py', '--clinSig', cls.clin_sig, '--ignore',
                cls.ignore, '--out', cls.out, '-e', cls.efo_map_file, '-g', cls.snp_2_gene_file,
                '--ot-schema', cls.ot_schema_path]
        cls.argparser = utilities.ArgParser(argv)

    def test_clin_sig(self):
        self.assertEqual(self.argparser.clinical_significance, self.clin_sig)

    def test_ignore(self):
        self.assertEqual(self.argparser.ignore_terms_file, self.ignore)

    def test_out(self):
        self.assertEqual(self.argparser.out, self.out)

    def test_efo_map_file(self):
        self.assertEqual(self.argparser.efo_mapping_file, self.efo_map_file)

    def test_snp_2_gene_file(self):
        self.assertEqual(self.argparser.snp_2_gene_file, self.snp_2_gene_file)

    def test_ot_schema_path(self):
        self.assertEqual(self.argparser.ot_schema, self.ot_schema_path)


class CheckDirExistsCreateTest(unittest.TestCase):
    def test_create(self):
        directory = "./test_tmp"
        utilities.check_dir_exists_create(directory)
        self.assertTrue(os.path.exists(directory))
        os.rmdir(directory)
