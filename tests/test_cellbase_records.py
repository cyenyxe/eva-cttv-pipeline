import unittest

from eva_cttv_pipeline import cellbase_records, config


class CellbaseRecordsTest(unittest.TestCase):
    def test_get_curr_result_list(self):
        curr_result_list = cellbase_records.get_curr_result_list(0)
        self.assertEqual(len(curr_result_list), config.BATCH_SIZE)

        curr_response = cellbase_records.get_curr_response(0)
        curr_result_list = cellbase_records.get_curr_result_list(curr_response['numTotalResults'])
        self.assertEqual(len(curr_result_list), 0)

        curr_result_list = cellbase_records.get_curr_result_list(999999)
        self.assertEqual(len(curr_result_list), 0)

        curr_response = cellbase_records.get_curr_response(0)
        len_to_expect = 20
        curr_result_list = cellbase_records.get_curr_result_list(curr_response['numTotalResults'] - len_to_expect)
        self.assertEqual(len(curr_result_list), len_to_expect)

    # def test_get_curr_result_lists(self):
    #     curr_response = clinvar_to_evidence_strings.get_curr_response(0)
    #     num_total_results = curr_response['numTotalResults']
    #     list_counter = 0
    #     for list in clinvar_to_evidence_strings.get_curr_result_lists():
    #         list_counter += 1
    #     pred_num_lists = math.ceil(num_total_results / config.BATCH_SIZE)
    #     self.assertEqual(pred_num_lists, list_counter)