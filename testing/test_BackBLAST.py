import unittest

from ..back_lib import run_blastp, create_hit_dict, HighScoringSegmentPair, filter_hsp, get_proteome_dict
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

		hsp_list = [HighScoringSegmentPair('ACY320321', 'ACY320320', 24, 1e-40, 60, 400),
		            HighScoringSegmentPair('ACY320322', 'ACY320320', 26, 1e-40, 60, 400),
		            HighScoringSegmentPair('ACY320323', 'ACY320320', 25, 1e-40, 60, 400),
		            HighScoringSegmentPair('ACY320324', 'ACY320320', 26, 1e-24, 60, 400),
		            HighScoringSegmentPair('ACY320325', 'ACY320320', 26, 1e-26, 60, 400),
		            HighScoringSegmentPair('ACY320326', 'ACY320320', 26, 1e-25, 60, 400),
		            HighScoringSegmentPair('ACY320326', 'ACY320320', 26, 0, 60, 400),
		            HighScoringSegmentPair('ACY320327', 'ACY320320', 26, 1e-40, 49, 400),
		            HighScoringSegmentPair('ACY320328', 'ACY320320', 26, 1e-40, 50, 400),
		            HighScoringSegmentPair('ACY320329', 'ACY320320', 26, 1e-40, 51, 400),
		            HighScoringSegmentPair('ACY320330', 'ACY320320', 26, 1e-40, 60, 299),
		            HighScoringSegmentPair('ACY320331', 'ACY320320', 26, 1e-40, 60, 300),
		            HighScoringSegmentPair('ACY320332', 'ACY320320', 26, 1e-40, 60, 301)]

		fail_accessions = ['ACY320321', 'ACY320324', 'ACY320327', 'ACY320330']

		for hsp in hsp_list:
			remove = filter_hsp(hsp, bitscore_cutoff=300)
			if remove:
				self.assertIn(hsp.query_seq_id, fail_accessions)
			else:
				self.assertNotIn(hsp.query_seq_id, fail_accessions)

	def test_get_proteome_dict(self):
		test_fasta_path = "./testing/test_files/test_databases/CP001220.1.faa"
		with open(test_fasta_path, 'r') as test_fasta:
			test_accession = test_fasta.readline().rsplit()[0].replace('>', '')
			test_line = test_fasta.readline()

			proteome_dict = get_proteome_dict(test_fasta_path)
			returned_test_fasta = proteome_dict[test_accession]
			test_fasta_line = returned_test_fasta.splitlines()[1]

			self.assertIn(test_fasta_line, test_line)

if __name__ == '__main__':
	unittest.main()
