import unittest

from ..back_lib import run_blastp, create_hit_dict
import csv


class MyTestCase(unittest.TestCase):
	def test_BLAST(self):
		with open("./testing/test_files/expected_outputs/forward_blast_one.txt", 'r') as blast_test_data:
			reference_blast_results = blast_test_data.read()
			reference_blast_results = reference_blast_results.replace(' ', '')
			new_blast_results = run_blastp("./testing/test_files/test_queries/CommamonasCluster.faa",
			                               "./testing/test_files/test_databases/CP001220.1.faa", "1e-25", 1)

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

			pass_test = True

			if first_hit.query_seq_id != test_query_seq_id:
				pass_test = False
			elif first_hit.sub_seq_id != test_sub_seq_id:
				pass_test = False
			elif first_hit.percent_identity != test_percent_identity:
				pass_test = False
			elif first_hit.evalue != test_evalue:
				pass_test = False
			elif first_hit.coverage != test_coverage:
				pass_test = False
			elif first_hit.bitscore != test_bitscore:
				pass_test = False

			self.assertTrue(pass_test)

if __name__ == '__main__':
	unittest.main()
