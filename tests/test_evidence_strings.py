import unittest

import eva_cttv_pipeline.evidence_strings as ES


class CTTVGeneticsEvidenceStringTest(unittest.TestCase):
    def setUp(self):
        self.test_ges = ES.CTTVGeneticsEvidenceString()

    def test_set_db_xref_url(self):
        url = "http://identifiers.org/clinvar.record/RCV000128628"
        self.test_ges.set_db_xref_url(url)
        self.assertEqual(self.test_ges['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'], url)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'], url)

    def test_set_url(self):
        url = "http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628"
        self.test_ges.set_url(url)
        self.assertEqual(self.test_ges['evidence']['gene2variant']['urls'][0]['url'], url)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['urls'][0]['url'], url)


# def test_dict_setter(setter_function, attrs_to_check, test_values):
#     for value in test_values:
#         setter_function(value)
#         for attr in attrs_to_check:
#             unittest.TestCase.assertEquals(value, attr)
#
#
# class CTTVGeneticsEvidenceStringTest(unittest.TestCase):
#     def setUp(self):
#         self.gen_ev_string = ES.CTTVGeneticsEvidenceString()
#
#     def test_set_db_xref_url(self):
#         test_urls = ["http://purl.obolibrary.org/obo/GO_0044691", "http://purl.obolibrary.org/obo/HP_0000113", "http://www.orpha.net/ORDO/Orphanet_71517"]
#         ev_urls = [self.gen_ev_string['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'], self.gen_ev_string['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url']]
#         test_dict_setter(self.gen_ev_string.set_db_xref_url, ev_urls, test_urls)
#         #
#         #
#         # for url in test_urls:
#         #     self.gen_ev_string.set_db_xref_url(url)
#         #     ev_url1 = self.gen_ev_string['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url']
#         #     self.assertEqual(url, ev_url1)
#         #     ev_url2 = self.gen_ev_string['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url']
#         #     self.assertEqual(url, ev_url2)
