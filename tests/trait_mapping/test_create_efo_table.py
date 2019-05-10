#!/usr/bin/env python3

import os
import sys
import tempfile
import unittest

sys.path.append('../../../')
from bin.trait_mapping.create_efo_table import *


class TestCreateEFOTableComponents(unittest.TestCase):

    def setUp(self):
        self.input_file = tempfile.NamedTemporaryFile(delete=False)
        self.output_file = tempfile.NamedTemporaryFile(delete=False)
        with open(self.input_file.name, 'w') as infile:
            infile.write('http://purl.obolibrary.org/obo/MONDO_0014263\n'
                         'http://www.orpha.net/ORDO/Orphanet_505216\n')

    def tearDown(self):
        os.remove(self.input_file.name)
        os.remove(self.output_file.name)

    def test_get_ols_details(self):
        """Check that the output of the script is adequate in terms of format; also check two example traits."""
        create_efo_table(self.input_file.name, self.output_file.name)
        with open(self.output_file.name) as outfile:
            out_lines = outfile.read().splitlines()
        self.assertEqual(len(out_lines), 2, 'Should get two output lines for two input lines')
        trait_1, trait_2 = [out_line.split('\t') for out_line in out_lines]
        for trait in [trait_1, trait_2]:
            self.assertEqual(len(trait), 22, 'Each output line must contain correct number of columns')

        # Trait 1
        self.assertEqual(trait_1[0], 'Verheij syndrome')  # disease
        self.assertEqual(trait_1[1], 'inherited genetic disease')  # child of
        self.assertEqual(trait_1[12], 'OMIM:615583')  # OMIM cross-ref
        self.assertEqual(trait_1[15], 'UMLS:C3810023')  # ULMS cross-ref
        self.assertEqual(trait_1[21], 'MONDO_0014263')  # MONDO ref

        # Trait 2
        self.assertEqual(trait_2[0], '3-methylglutaconic aciduria type 9')  # disease
        self.assertEqual(trait_2[1], 'disease')  # child of
        self.assertEqual(trait_2[3], '3-methylglutaconic aciduria-epilepsy-spasticity-severe intellectual '
                                     'disability syndrome || MGA9')  # synonyms
        self.assertEqual(trait_2[18], 'Orphanet:505216')  # Orphanet ref
