#!/usr/bin/env python3

import re
import os.path
import sys

ltrans = {1: 'Class', 2: 'Mechanism', 3: 'Group', 4: 'Model'}

reads_mapped = 0
total_reads = 0


def load_model_members(infile):
	return_values = {}
	with open(infile, 'r') as f:
		data = f.read().split('\n')
		for line in data[1:]:
			if not line:
				continue
			entries = line.split('\t')
			return_values.setdefault(entries[0], [x[1:] if x[0] == '@' else x for x in entries[1].split(',')])
	return {y: k for k, v in return_values.items() for y in v}


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


def load_model_annotations(infile):
	return_values = {}
	with open(infile, 'r') as f:
		data = f.read().split('\n')
		for line in data[1:]:
			if not line:
				continue
			entries = line.split('\t')
			return_values.setdefault(entries[1], entries[4:] + [entries[2]])
	return return_values


def find_substring(sub, s_array):
	found = ""
	for d in s_array:
		if len(sub) > len(d):
			continue
		if sub in d:
			return d
		for i in range(len(d) - len(sub) + 1):
			if d[i:i+len(sub)] == sub:
				return d
	if not found:
		sys.stderr.write("No match for {}\n".format(sub))
		raise



def load_tblout(infile):
	ret = {}
	with open(infile, 'r') as f:
		data = f.read().split('\n')
		for line in data:
			if line.startswith('#') or not line:
				continue
			entry = line.split()
			contig = entry[0]
			if int(entry[6]) < int(entry[7]):
				ali = (int(entry[6]), int(entry[7]))
			else:
				ali = (int(entry[7]), int(entry[6]))
			ret.setdefault(contig, []).append((ali[0], ali[1], entry[2]))
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
					return temp[2], temp[0], int(temp[3]), temp[5]  # RefName, readname, 1-start, CIGAR

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
	counts = [{}, {}, {}, {}]
	output_annots = []

	D = [load_tblout(x) for x in sys.argv[1:4]]
	A = load_annotations(sys.argv[4])
	MA = load_model_annotations(sys.argv[5])
	MM = load_model_members(sys.argv[6])
	param_string = sys.argv[1].split('/')[-2].replace('_', ',')
	total = sys.argv[7]
	for refname, simulated, start, cigar in SamParser('-'):
		if simulated[0:2] == 'gi':
			continue
		# RefName, 1-start, CIGAR, RefLen, ReadSeq
		stop = start + parse_cigar(cigar) - 1
		for d in D:
			if refname not in d:
				continue
			model_hits = set()
			for triplets in d[refname]:
				if max(start, triplets[0]) <= min(stop, triplets[1]):
					if MA[triplets[2]] != 'NA':
						model_hits.add(triplets[2])
			header = '-'.join(simulated.split('-')[:-2])
			annots = A[header]
			output_annots = [x for x in annots]
			if header not in MM:
				true_model = find_substring(header.replace('|RequiresSNPConfirmation', '').replace('CARD|', ''), MM.keys())
				MM.setdefault(header, MM[true_model])
			target_model = MM[header]
			output_annots.append(target_model)
			if model_hits:
				for target in model_hits:
					target_annots = MA[target][:-1]
					for a in range(len(annots)):
						if annots[a] in target_annots[a].split('|'):
							counts[a].setdefault(simulated, 0)
							if counts[a][simulated] < 2:
								counts[a][simulated] += 1
						if target == target_model:
							counts[3].setdefault(simulated, 0)
							if counts[3][simulated] < 2:
								counts[3][simulated] += 1
	for level in range(len(output_annots)):
		sys.stdout.write('mmarc_assembly,{},{},{},{},{}\n'.format(
			param_string,
			level,
			output_annots[level],
			str(sum([y for y in counts[level].values()])),
			total
		))
