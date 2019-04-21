#!/usr/bin/env python3

import sys


def count_sam_headers_simulated(infile):
	with open(infile, 'r') as f:
		counter = 0
		line = f.readline()
		while line:
			if (not line[0] == '@') and (not line[0:2] == 'gi'):
				counter += 1
			line = f.readline()
	sys.stdout.write(str(counter))


if __name__ == '__main__':
	count_sam_headers_simulated(sys.argv[1])
