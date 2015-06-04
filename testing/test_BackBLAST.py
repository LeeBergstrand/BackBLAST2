import unittest

from ..back_lib import run_blastp


class MyTestCase(unittest.TestCase):
	def test_BLAST(self):
		with open("./testing/test_files/expected_outputs/forward_blast_one.txt",'r') as blast_test_data:
			reference_blast_results = blast_test_data.read()
			reference_blast_results = reference_blast_results.replace(' ', '')
			new_blast_results = run_blastp("./testing/test_files/test_queries/CommamonasCluster.faa", "./testing/test_files/test_databases/CP001220.1.faa", "1e-25", 1)

			self.assertEqual(new_blast_results, reference_blast_results)

if __name__ == '__main__':
	unittest.main()
