#!/usr/bin/env python3

import sys
import re


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


def count_sam_headers_simulated(infiles, A, MA, MM, total, params):
	counts = [{}, {}, {}, {}]
	output_annots = []
	for file_idx in range(len(infiles)):
		with open(infiles[file_idx], 'r') as f:
			# class, mech, group
			line = f.readline()
			while line:
				while line[0] == '#':
					line = f.readline()
					if not line:
						break
				if not line:
					break
				if not line[0:2] == 'gi':
					simulated = line.split()[0].split('/')[0]
					header = '-'.join(simulated.split('-')[:-2])
					target = line.split()[2]
					annots = A[header]
					output_annots = [x for x in annots]
					target_annots = MA[target][:-1]
					if header not in MM:
						true_model = find_substring(header.replace('|RequiresSNPConfirmation', '').replace('CARD|', ''), MM.keys())
						MM.setdefault(header, MM[true_model])
					target_model = MM[header]
					output_annots.append(target_model)
					for a in range(len(annots)):
						if annots[a] in target_annots[a].split('|'):
							counts[a].setdefault(simulated, 0)
							if counts[a][simulated] < 2:
								counts[a][simulated] += 1
						# else:
						# 	sys.stdout.write('{}\t{}\n'.format(annots[a], target_annots[a]))
						if target == target_model:
							counts[3].setdefault(simulated, 0)
							if counts[3][simulated] < 2:
								counts[3][simulated] += 1
				line = f.readline()
	for level in range(len(output_annots)):
		sys.stdout.write('mmarc_hts,{},{},{},{},{}\n'.format(
			params,
			level,
			output_annots[level],
			str(sum([y for y in counts[level].values()])),
			total
		))


if __name__ == '__main__':
	# hts_tbloutI, hts_tbloutII, hts_tbloutIII, megares_annotations, mmarc_annotations, total
	param_string = sys.argv[1].split('/')[-2].replace('_', ',')
	A = load_annotations(sys.argv[4])
	MA = load_model_annotations(sys.argv[5])
	MM = load_model_members(sys.argv[6])
	count_sam_headers_simulated([sys.argv[1], sys.argv[2], sys.argv[3]], A, MA, MM, int(sys.argv[7]), param_string)
