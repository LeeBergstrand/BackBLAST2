import csv
import subprocess
import sys
from hashlib import md5
from Bio import SeqIO

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
def run_blastp(query_file, database_file, e_value_cutoff, processes):
	"""
	Runs BLASTp...
	:param query_file: The amino acid query FASTA file.
	:param database_file: The amino acid BLAST database location (amino acid FASTA file at this location)
	:param e_value_cutoff: The e-value cutoff for BLASTp.
	:param processes: Number of processes for BLASTp to use.
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

# 3: Filters HSPs by Percent Identity...
def filter_blast_csv(BLASTOut):
	minIdent = 25

	BLASTCSVOut = BLASTOut.splitlines(True)  # Converts raw BLAST csv output into list of csv rows.
	BLASTreader = csv.reader(BLASTCSVOut)  # Reads BLAST csv rows as a csv.

	BLASTCSVOutFiltred = []  # Note should simply delete unwanted HSPs from current list rather than making new list.
	# Rather than making a new one.
	for HSP in BLASTreader:
		if HSP[2] >= minIdent:  # Filters by minimum identity.
			# Converts each HSP parameter that should be a number to a number.
			HSP[2] = float(HSP[2])
			HSP[3] = float(HSP[3])
			HSP[4] = float(HSP[4])
			HSP[5] = float(HSP[5])
			BLASTCSVOutFiltred.append(HSP)  # Appends to output array.

	return BLASTCSVOutFiltred


# -------------------------------------------------------------------------------------------------
# 5: Creates a python dictionary (hash table) that contains the the FASTA for each protein in the proteome.
def createProteomeHash(ProteomeFile):
	ProteomeHash = dict()
	try:
		handle = open(ProteomeFile, "rU")
		proteome = SeqIO.parse(handle, "fasta")
		for record in proteome:
			ProteomeHash.update({record.id: record.format("fasta")})
		handle.close()
	except IOError:
		print("Failed to open " + ProteomeFile)
		sys.exit(1)

	return ProteomeHash
