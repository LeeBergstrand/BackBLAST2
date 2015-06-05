import unittest

from ..back_lib import run_blastp, create_hit_dict, HighScoringSegmentPair, filter_hsp
import csv


class MyTestCase(unittest.TestCase):
	def test_BLAST(self):
		with open("./testing/test_files/expected_outputs/forward_blast_one.txt", 'r') as blast_test_data:
			reference_blast_results = blast_test_data.read()
			reference_blast_results = reference_blast_results.replace(' ', '')
			new_blast_results = run_blastp("./testing/test_files/test_queries/CommamonasCluster.faa",
			                               "./testing/test_files/test_databases/CP001220.1.faa", 1, "1e-25")

			self.assertEqual(new_blast_results, reference_blast_results)

	def test_hit_loading(self):
		with open("./testing/test_files/expected_outputs/forward_blast_one_hit.txt", 'r') as blast_test_data:
			raw_hit_data = blast_test_data.read()
			reference_blast_results = raw_hit_data.replace(' ', '')

			blast_data_rows = reference_blast_results.splitlines(True)
			blast_csv = csv.reader(blast_data_rows)

			row = next(blast_csv)
			test_query_seq_id = str(row[0])
			test_sub_seq_id = str(row[1])
			test_percent_identity = float(row[2])
			test_evalue = float(row[3])
			test_coverage = float(row[4])
			test_bitscore = float(row[5])

			hit_dict = create_hit_dict(reference_blast_results)
			first_hit = hit_dict.popitem()[1]

			self.assertTrue(first_hit.query_seq_id == test_query_seq_id)
			self.assertTrue(first_hit.sub_seq_id == test_sub_seq_id)
			self.assertTrue(first_hit.percent_identity == test_percent_identity)
			self.assertTrue(first_hit.evalue == test_evalue)
			self.assertTrue(first_hit.coverage == test_coverage)
			self.assertTrue(first_hit.bitscore == test_bitscore)

	def test_hit_filtering(self):

		hsp_list = []
		hsp_list.append(HighScoringSegmentPair('ACY320321', 'ACY320320', 24, 1e-40, 60, 400))  # Fail (identity)
		hsp_list.append(HighScoringSegmentPair('ACY320322', 'ACY320320', 26, 1e-40, 60, 400))  # Pass
		hsp_list.append(HighScoringSegmentPair('ACY320323', 'ACY320320', 25, 1e-40, 60, 400))  # Pass
		hsp_list.append(HighScoringSegmentPair('ACY320324', 'ACY320320', 26, 1e-24, 60, 400))  # Fail (evalue)
		hsp_list.append(HighScoringSegmentPair('ACY320325', 'ACY320320', 26, 1e-26, 60, 400))  # Pass
		hsp_list.append(HighScoringSegmentPair('ACY320326', 'ACY320320', 26, 1e-25, 60, 400))  # Pass
		hsp_list.append(HighScoringSegmentPair('ACY320326', 'ACY320320', 26, 0, 60, 400))      # Pass
		hsp_list.append(HighScoringSegmentPair('ACY320327', 'ACY320320', 26, 1e-40, 49, 400))  # Fail (coverage)
		hsp_list.append(HighScoringSegmentPair('ACY320328', 'ACY320320', 26, 1e-40, 50, 400))  # Pass
		hsp_list.append(HighScoringSegmentPair('ACY320329', 'ACY320320', 26, 1e-40, 51, 400))  # Pass
		hsp_list.append(HighScoringSegmentPair('ACY320330', 'ACY320320', 26, 1e-40, 60, 299))  # Fail (bitscore)
		hsp_list.append(HighScoringSegmentPair('ACY320331', 'ACY320320', 26, 1e-40, 60, 300))  # Pass
		hsp_list.append(HighScoringSegmentPair('ACY320332', 'ACY320320', 26, 1e-40, 60, 301))  # Pass

		fail_accessions = ['ACY320321', 'ACY320324', 'ACY320327', 'ACY320330']

		for hsp in hsp_list:
			remove = filter_hsp(hsp, bitscore_cutoff=300)
			if remove:
				self.assertIn(hsp.query_seq_id, fail_accessions)
			else:
				self.assertNotIn(hsp.query_seq_id, fail_accessions)


if __name__ == '__main__':
	unittest.main()
