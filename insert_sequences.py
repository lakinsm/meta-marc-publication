#!/usr/bin/env python3

import sys
import random
import math
random.seed(154)

NUCLEOTIDES = ('A', 'C', 'G', 'T')


def fastaParse(infile):
	with open(infile, 'r') as fasta_file:
		# Skip whitespace
		while True:
			line = fasta_file.readline()
			if line is "":
				return  # Empty file or premature end of file?
			if line[0] is ">":
				break
		while True:
			if line[0] is not ">":
				raise ValueError("Records in FASTA should begin with '>'")
			header = line[1:].rstrip()
			all_lines = []
			line = fasta_file.readline()
			while True:
				if not line:
					break
				if line[0] is ">":
					break
				all_lines.append(line.rstrip())
				line = fasta_file.readline()
			yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
			if not line:
				return # Stop Iteration


def selectGeneAndInsert(database, gene_number, copy_number, insertion_size, contiguous, out_dir):
	random_gene_idxs = random.sample(range(len(database)), int(gene_number))
	for idx in random_gene_idxs:
		header, seq = database[idx]
		bp_to_replace = math.ceil(float(insertion_size) * len(seq))
		bp_start_idx = random.choice(range(len(seq) - bp_to_replace - 1))
		out_seq = ""
		if contiguous == "true":
			out_seq += seq[:bp_start_idx]
			for i in range(bp_start_idx, bp_start_idx + bp_to_replace):
				out_seq += random.choice(NUCLEOTIDES)
			out_seq += seq[bp_start_idx + bp_to_replace:]
		else:
			random_nucleotide_idxs = set(random.sample(range(len(seq)), bp_to_replace))
			for i in range(len(seq)):
				if i in random_nucleotide_idxs:
					out_seq += random.choice(NUCLEOTIDES)
				else:
					out_seq += seq[i]
		with open(out_dir + '/' + '_'.join(
				[str(x) for x in (idx, insertion_size, copy_number, contiguous)]
			) + ".fasta", 'w') as out:
			for copy in range(int(copy_number)):
				out.write('>{}-copy{}\n{}\n'.format(
					header, copy, out_seq
				))


if __name__ == '__main__':
	# 1: MEGARes database filepath
	# 2: Number of genes to randomly sample
	# 3: Number of copies to make of modified gene
	# 4: Proportion of gene to replace (Y)
	# 5: Contiguous mutations [true, false]
	# 6: Path of output directory
	D = [(k, v) for k, v in fastaParse(sys.argv[1])]
	selectGeneAndInsert(D, sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
