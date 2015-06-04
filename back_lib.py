import csv
import subprocess
import sys
from Bio import SeqIO


# Functions:
# =================================================================================================
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
# 3: Filters HSPs by Percent Identity...
def filterBLASTCSV(BLASTOut):
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
