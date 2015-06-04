import unittest

from ..back_lib import run_blastp


class MyTestCase(unittest.TestCase):
	def test_BLAST(self):
		results = run_blastp("./testing/test_files/test_queries/CommamonasCluster.faa", "./testing/test_files/test_databases/CP001220.1.faa", "1e-25", 1)
		self.assertTrue(results)

if __name__ == '__main__':
	unittest.main()
