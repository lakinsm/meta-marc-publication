#!/usr/bin/env python3

import sys
import os.path
import re

read_mapped = 0
total_reads = 0


def load_annotations(infile):
	ret = {}
	with open(infile, 'r') as f:
		data = f.read().split('\n')
		for line in data:
			if not line:
				continue
			entries = line.split(',')
			ret.setdefault(entries[0], entries[1:])
	return ret


def load_resfams_metadata(infile):
	ret = {}
	with open(infile, 'r') as f:
		data = f.read().split('\n')[1:]
		for line in data:
			if not line:
				continue
			entry = line.split('\t')
			ret[entry[0]] = entry[1]
	return ret


def load_domains(infile):
	ret = {}
	with open(infile, 'r') as f:
		data = f.read().split('\n')
		for line in data:
			if line.startswith('#') or not line:
				continue
			entry = line.split()
			contig = entry[-3].replace('>', '')
			locs = [int(x) for x in entry[0].split('|')[-2:]]
			ret.setdefault(contig, []).append((int(locs[0]), int(locs[1]), entry[3]))
	return ret


def parse_cigar(s):
	length = 0
	ret = re.findall(r'(\d+)([A-Z=]{1})', s)
	universe = {'X', 'P', 'I', 'N', 'D', '=', 'M'}
	for occ, op in ret:
		if op in universe:
			length += int(occ)
	return length


class SamParser:
	"""This object takes as input a SAM file path and constructs an iterable that outputs
    hash-mapping of header to sequence information.  Only one line will be held in memory at a time using this method.
    """

	def __init__(self, filepath):
		"""
        constructor
        @param filepath: filepath to the input raw SAM file.
        """
		if os.path.exists(filepath):  # if file is a file, read from the file
			self.sam_file = str(filepath)
			self.stdin = False
		elif not sys.stdin.isatty():  # else read from standard in
			self.stdin = True
		else:
			raise ValueError("Parameter filepath must be a SAM file")
		self.current_line = None
		self.reads_mapping = 0
		self.reads_total = 0
		self.header_lens = {}

	def __iter__(self):
		return self

	@property
	def _iterate(self):
		# Skip all leading whitespace
		while True:
			if self.stdin:
				sam_line = sys.stdin.readline()  # read from stdin
			else:
				sam_line = self.sam_file.readline()  # read from file
			if not sam_line:
				return  # End of file
			if sam_line[0] != '@':  # these lines are the actual reads
				self.reads_total += 1
				temp = sam_line.split()
				if (int(temp[1]) & 4) == 0:
					self.reads_mapping += 1
					return temp[2], temp[0], int(temp[3]), temp[5]  # RefName, header, 1-start, CIGAR
		self.sam_file.close()  # catch all in case this line is reached
		assert False, "Should not reach this line"

	def __next__(self):
		if not self.stdin and type(self.sam_file) is str:  # only open file here if sam_file is a str and not fileIO
			self.sam_file = open(self.sam_file, "r")
		value = self._iterate
		if not value:  # close file on EOF
			if not self.stdin:
				self.sam_file.close()
			global reads_mapped
			global total_reads
			reads_mapped = self.reads_mapping
			total_reads = self.reads_total
			raise StopIteration()
		else:
			return value


if __name__ == '__main__':
	counts = {}
	output_annot = ""

	D = load_domains(sys.argv[1])
	A = load_annotations(sys.argv[2])
	R = load_resfams_metadata(sys.argv[3])
	total = sys.argv[4]
	param_string = sys.argv[1].split('/')[-2].replace('_', ',')

	for refname, simulated, start, cigar in SamParser('-'):
		if simulated[0:2] == 'gi':
			continue
		# RefName, 1-start, CIGAR, RefLen, ReadSeq
		stop = start + parse_cigar(cigar) - 1
		header = '-'.join(simulated.split('-')[:-2])
		class_annot = A[header][0]
		output_annot = class_annot
		if refname not in D:
			continue
		model_hits = set()
		for triplets in D[refname]:
			if max(start, triplets[0]) <= min(stop, triplets[1]):
				if R[triplets[2]] != 'NA':
					model_hits.add(triplets[2])
		if model_hits:
			for target in model_hits:
				class_target = R[target]
				if class_annot == class_target:
					counts.setdefault(simulated, 0)
					if counts[simulated] < 2:
						counts[simulated] += 1
	sys.stdout.write('resfams,{},{},{},{},{}\n'.format(
		param_string,
		0,
		output_annot,
		str(sum([y for y in counts.values()])),
		total
	))
