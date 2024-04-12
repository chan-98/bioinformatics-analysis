#!/usr/bin/env python3

"""
Chandini
12/04/2024

Script to replace IUPAC alleles in consensus with max-depth base from mpileup and create new consensus.
Usage:
        python3 generate_new_consensus.py consensus.fasta mpileup.txt consensus_replaced.fasta
"""

import os
import sys
from Bio import SeqIO

class Consensus:

	def __init__(self, record):
		self.record = record
		self.id = record.id
		self.seq = record.seq
	
	def get_iupac_positions(self):
		positions = []
		bases = ['a','t','g','c','n','A','T','G','C','N']
		for i,b in enumerate(self.seq):
			if b not in bases:
				positions.append(i)
		return positions
		
	def get_max_depth_base_at(self,pos,mpileup_file):
		pos = {i+1:'n' for i in pos} # 1 based indexing
		iupac = { "K":["G", "T"], "R":["A", "G"], "M":["A", "C"], "S":["G", "C"], "Y":["T","C"], "W":["A","T"] }

		with open(mpileup_file) as f:
			for line in f:
				row = line.rstrip().split('\t')
				if int(row[1]) in pos.keys():
					bases = [row[2], 'A', 'T', 'G', 'C', 'N']
					A = row[4].count('A')+row[4].count('a')
					T = row[4].count('T')+row[4].count('t')
					G = row[4].count('G')+row[4].count('g')
					C = row[4].count('C')+row[4].count('c')
					ref = row[4].count('.')+row[4].count(',')
					base_count = [ref,A,T,G,C]
					ix = base_count.index(max(base_count))
					pos[int(row[1])] = bases[ix]
		return list(pos.values())


def write_new_consensus(ref, mpileup, new_consensus):

	headers = []
	sequences = []

	for record in SeqIO.parse(ref,"fasta"):

		print(f"Writing new consensus sequence for the record {record.id}...")
		c = Consensus(record)
		positions = c.get_iupac_positions()
		bases = c.get_max_depth_base_at(positions, mpileup)
			
		new_seq = ""
		sample_name = os.path.basename(new_consensus).split(".")[0]
			
		for i,b in enumerate(record.seq):
			if i in positions:
				ix = positions.index(i)
				new_base = bases[ix]
				new_seq += new_base
			else:
				new_seq += b
	
		headers.append(">"+record.id+"|"+sample_name)
		sequences.append(new_seq)

	print("Writing Output File")
	with open(new_consensus,"w") as f:
		for i in range(0, len(headers)):
			f.write(headers[i]+"\n")
			f.write(sequences[i]+"\n")


if __name__ == "__main__":

	if len(sys.argv) != 4:
		print("USAGE: python3", sys.argv[0], "[ref.fasta] [mpileup.txt] [output.fasta]")
	else:
		reference = sys.argv[1]
		mpileup = sys.argv[2]
		output = sys.argv[3]
		write_new_consensus(reference, mpileup, output)
		print("FINISHED.")
