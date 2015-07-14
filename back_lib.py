#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2015)

Description: Python library for the BackBLAST reciprocal BLAST program.

Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
              - All operations are to be done with protein sequences.
"""

import csv
import subprocess
import sys
from hashlib import md5
from Bio import SeqIO

# ==============================
# Python 2/3 Compatibility Code:
# ==============================

# Creates custom functions for Python2/3 dictionary iteration.
try:
	dict.iteritems
except AttributeError:
	# Python 3
	def itervalues(d):
		return iter(d.values())

	def iteritems(d):
		return iter(d.items())
else:
	# Python 2
	def itervalues(d):
		return d.itervalues()

	def iteritems(d):
		return d.iteritems()

# ==========
# Classes:
# ==========

# -------------------------------------------------------------------------------------------------
class HighScoringSegmentPair(object):
	"""
	Definition of a HighScoringSegmentPair.
	:var query_seq_id: The query protein's sequence ID
	:var sub_seq_id: The query protein's sequence ID
	:var percent_identity: The percentage of identical amino acids between the two proteins.
	:var evalue:    The expect value describes the number of hits one can "expect" to see by chance when searching
					a database of a particular size.
	:var coverage: The amount of aligned sequence overlap between the two proteins.
	:var bitscore: A score of the quality of alignment between the two proteins which takes into account database size.
	"""

	def __init__(self, query_seq_id, sub_seq_id, percent_identity, evalue, coverage, bitscore):
		self.query_seq_id = str(query_seq_id)
		self.sub_seq_id = str(sub_seq_id)
		self.percent_identity = float(percent_identity)
		self.evalue = float(evalue)
		self.coverage = float(coverage)
		self.bitscore = float(bitscore)

	def get_md5(self):
		hash_string = "".join([str(x) for x in self.__dict__.values()])  # Join all attributes into a single string.
		hash_md5 = md5(hash_string.encode('utf-8')).hexdigest()  # Create md5 hash.
		return hash_md5


# ==========
# Functions:
# ==========

# -------------------------------------------------------------------------------------------------
def run_blastp(query_file, database_file, processes, e_value_cutoff=1e-25):
	"""
	Runs BLASTp...
	:param query_file: The amino acid query FASTA file.
	:param database_file: The amino acid BLAST database location (amino acid FASTA file at this location)
	:param processes: Number of processes for BLASTp to use.
	:param e_value_cutoff: The e-value cutoff for BLASTp (Default = 1e-25).
	:return:    A csv formatted BLASTp output (query_sequence_id, subject_sequence_id, percent_identity, e-value,
				query coverage, bitscore)
	"""
	blast_out = subprocess.check_output(
		["blastp", "-db", database_file, "-query", query_file, "-evalue", str(e_value_cutoff), "-num_threads",
		 str(processes), "-outfmt", "10 qseqid sseqid pident evalue qcovhsp bitscore"])

	# Decodes BLASTp output to UTF-8 (In Py3 check_output returns raw bytes)
	blast_out = blast_out.decode().replace(' ', '')
	return blast_out


# -------------------------------------------------------------------------------------------------
def create_hit_dict(blast_csv_in):
	"""
	Creates a dictionary of HSP objects.
	:param blast_csv_in: CSV formatted multi-line string BLASTp results from the runBLAST method.
	:return: Dictionary of HSP objects with a md5 hash of each object's attributes as key.
	"""
	hit_dict = {}

	blast_csv_in = blast_csv_in.splitlines(True)  # Converts raw BLAST csv output into list of csv rows.
	blast_reader = csv.reader(blast_csv_in)  # Reads BLAST csv rows as a csv.

	for row in blast_reader:
		query_seq_id = row[0]
		sub_seq_id = row[1]
		percent_identity = row[2]
		evalue = row[3]
		coverage = row[4]
		bitscore = row[5]

		new_hsp = HighScoringSegmentPair(query_seq_id, sub_seq_id, percent_identity, evalue, coverage, bitscore)
		hit_dict[new_hsp.get_md5()] = new_hsp

	return hit_dict


# -------------------------------------------------------------------------------------------------
def filter_hsp_dict(raw_hit_dict, percent_identity_cutoff=25, evalue_cutoff=1e-25, coverage_cutoff=50,
                    bitscore_cutoff=0):
	"""
	Filters HSP dict by removing HSPs with percent identity, e-value, coverage, or bitscore below a cutoff value.
	:param raw_hit_dict: Input dictionary of HSP objects.
	:param percent_identity_cutoff: The percent identity cutoff (Default = 25%).
	:param evalue_cutoff: The e-value cutoff (Default = 1e-25).
	:param coverage_cutoff: The percent alignment coverage cutoff (Default = 50%).
	:param bitscore_cutoff: The bitscore cutoff (Default = 0).
	:return: A filtered dict of hsp objects.
	"""

	filtered_hsp_dict = {}

	for key, hsp in iteritems(raw_hit_dict):
		if not filter_hsp(hsp, percent_identity_cutoff, evalue_cutoff, coverage_cutoff, bitscore_cutoff):
			filtered_hsp_dict[key] = hsp

	return filtered_hsp_dict


# -------------------------------------------------------------------------------------------------
def filter_hsp(hsp, percent_identity_cutoff=25, evalue_cutoff=1e-25, coverage_cutoff=50, bitscore_cutoff=0):
	"""
	Filters by Percent Identity, e-value, coverage, and bitscore.
	:param hsp: Input dictionary of HSP objects.
	:param percent_identity_cutoff: The percent identity cutoff (Default = 25%).
	:param evalue_cutoff: The e-value cutoff (Default = 1e-25).
	:param coverage_cutoff: The percent alignment coverage cutoff (Default = 50%).
	:param bitscore_cutoff: The bitscore cutoff (Default = 0).
	:return: A boolean which indicates whether an HSP should be rejected.
	"""

	reject = False
	percent_identity = hsp.percent_identity
	evalue = hsp.evalue
	coverage = hsp.coverage
	bitscore = hsp.bitscore

	if float(percent_identity) < float(percent_identity_cutoff):
		reject = True
	elif float(evalue) > float(evalue_cutoff):
		reject = True
	elif float(coverage) < float(coverage_cutoff):
		reject = True
	elif float(bitscore) < float(bitscore_cutoff):
		reject = True

	return reject


# -------------------------------------------------------------------------------------------------
def get_proteome_dict(proteome_file):
	"""
	Creates a python dictionary (hash table) that contains the the FASTA for each protein in the proteome.
	:param proteome_file: An input amino acid FASTA file containing the CDS of an organism.
	:return:    A dictionary of amino acid FASTA sequences containing the CDS of an organism.
				The sequence ID of each CDS is used as a key.
	"""
	proteome_hash = dict()
	try:
		handle = open(proteome_file, "rU")
		proteome = SeqIO.parse(handle, "fasta")
		for record in proteome:
			proteome_hash.update({record.id: record.format("fasta")})
		handle.close()
	except IOError:
		print("Failed to open " + proteome_file)
		sys.exit(1)

	return proteome_hash
