#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2015)

Description: A program that takes query proteins and uses BLASTp to search for highly similar
             proteins within a local BLAST database derived from a target proteome. The program
             then uses BLASTp again to do a reverse search for found subject proteins in the query.
             It then filters hits down to Best Reciprocal Hits to confirm gene orthology.

Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
              - This script requires BLAST+ 2.2.9 or later.
              - All operations are to be done with protein sequences.
              - All query proteins should be from sequenced genomes in order to facilitate reciprocal BLASTp.
              - MakeAABlastDB must be used to create BLASTn databases for both query and subject proteomes.
              - BLASTp requires that the FASTA file the subject database remain in the same directory as the database.
"""

# Imports:
import argparse


def main(args):
	print("Hello World")


if __name__ == '__main__':
	# ------------------------------
	"""
	Command Line Interface Options:

	Type BackBLAST2.py -h for help.
	"""
	# ------------------------------

	parser = argparse.ArgumentParser()
	parser.add_argument('-q', '--query', metavar='QUERY', nargs=1, help='''
	The path to a FASTA file to be used as a query.''')

	parser.add_argument('-s', '--subject_database', metavar='DB', nargs=1, help='''
	The path to the subject BLAST database. ''')

	parser.add_argument('-r', '--query_database', metavar='DB', nargs=1, help='''
	The path to the subject BLAST database. ''')

	parser.add_argument('-o', '--outdir', metavar='OUTPATH', nargs=1, help='''
	An output directory for output files. If not specified, the current working directory is used.''')

	parser.add_argument('-j', '--jobs', metavar='JOBS', nargs=1, type=int, help='''
	The number of processes to use for multiple files (defaults to the number of processor cores).''')

	cli_args = parser.parse_args()

	# At minimum we require a query, query DB and subject DB to proceed.
	proceed = True

	if cli_args.query is None:
		print("Error: Missing query FASTA path...")
		proceed = False

	if cli_args.subject_database is None:
		print("Error: Missing subject BLAST database path...")
		proceed = False

	if cli_args.query_database is None:
		print("Error: Missing query BLAST database path...")
		proceed = False

	if proceed:
		main(cli_args)
	else:
		print("")
		parser.print_help()
