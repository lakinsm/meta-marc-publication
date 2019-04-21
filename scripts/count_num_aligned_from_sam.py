#!/usr/bin/env python3

import sys
import re


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


def count_sam_headers_simulated(infiles, A, total, params):
	counts = [{}, {}, {}]
	output_annots = []
	for file_idx in range(len(infiles)):
		with open(infiles[file_idx], 'r') as f:
			# class, mech, group
			line = f.readline()
			while line:
				if (not line[0] == '@') and (not line[0:2] == 'gi'):
					if int(line.split()[1]) & 4 == 0:
						simulated = line.split()[0]
						header = '-'.join(simulated.split('-')[:-2])
						target = line.split()[2]
						if target[0] == '@':
							target = target[1:]
						annots = A[header]
						output_annots = annots
						if target not in A:
							true_target = find_substring(target, A.keys())
							A.setdefault(target, A[true_target])
						target_annots = A[target]
						for a in range(len(annots)):
							if annots[a] == target_annots[a]:
								counts[a].setdefault(simulated, 0)
								if counts[a][simulated] < 2:
									counts[a][simulated] += 1
							# else:
							# 	sys.stdout.write('{}\t{}\n'.format(annots[a], target_annots[a]))
				line = f.readline()
	for level in range(len(output_annots)):
		sys.stdout.write('alignment,{},{},{},{},{}\n'.format(
			params,
			level,
			output_annots[level],
			str(sum([y for y in counts[level].values()])),
			total
		))


if __name__ == '__main__':
	# samI, samII, samIII, annotations, total
	param_string = sys.argv[1].split('/')[-2].replace('_', ',')
	A = load_annotations(sys.argv[4])
	count_sam_headers_simulated([sys.argv[1], sys.argv[2], sys.argv[3]], A, int(sys.argv[5]), param_string)
