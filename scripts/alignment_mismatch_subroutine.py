#!/usr/bin/env python3


import sys
import os.path
import re
import numpy as np


cov_vectors = {}
mismatches = {'Class': {}, 'Mechanism': {}, 'Group': {}, 'Gene': {}}
base_mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
level_map = {0: 'Class', 1: 'Mechanism', 2: 'Group'}

reads_mapped = 0
total_reads = 0


def load_annotations(infile):
    ret = {}
    with open(infile, 'r') as f:
        data = f.read().split('\n')[1:]
        for entry in data:
            if entry:
                line = entry.split(',')
                ret.setdefault(line[0], line[1:])
    return ret


def fasta_parse(infile):
    """ Parses a fasta file in chunks of 2 lines.
    :param infile: path to the input fasta file
    :return: generator of (header, sequence) fasta tuples
    """
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
                return  # Stop Iteration
        assert False, "Should not reach this line"


def parse_cigar_reference(s, refseq, ref_start, ref_name, readseq):
    ret = re.findall(r'(\d+)([A-Z=]{1})', s)
    if ref_name not in cov_vectors:
        cov_vectors[ref_name] = [np.zeros((4, len(refseq)), dtype=np.int32), 0, 0]
    cov_vectors[ref_name][2] += 1
    current_loc = ref_start
    for occ, op in ret:  # Number, operator
        if op in {'P', 'N'}:
            cov_vectors[ref_name][1] += 1
            current_loc += int(occ)
        elif op == 'I':
            cov_vectors[ref_name][1] += 1
        elif op == 'D':
            current_loc += int(occ)
        elif op in {'M', '=', 'X'}:
            for i in range(int(occ)):
                #print(ref_start, len(refseq), current_loc, occ, op, s)
                try:
                    cov_vectors[ref_name][0][base_mapping[readseq[i]]][current_loc] += 1
                    current_loc += 1
                except KeyError:
                    current_loc += 1
        elif op in {'S', 'H'}:
            pass




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
            if sam_line[0] != "@":  # these lines are the actual reads
                self.reads_total += 1
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rReads processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = sam_line.split()
                if (int(temp[1]) & 4) == 0:
                    self.reads_mapping += 1
                    return temp[2], int(temp[3])-1, temp[5], int(temp[8]), temp[9]  # RefName, 0-start, CIGAR, RefLen, ReadSeq
        self.sam_file.close()  # catch all in case this line is reached
        assert False, "Should not reach this line"

    def __next__(self):
        if not self.stdin and type(self.sam_file) is str:  # only open file here if sam_file is a str and not fileIO
            self.sam_file = open(self.sam_file, "r")
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.sam_file.close()
            sys.stdout.write("\n{:d} reads mapped out of {:d} total reads\n".format(self.reads_mapping, self.reads_total))
            sys.stdout.flush()
            global reads_mapped
            global total_reads
            reads_mapped = self.reads_mapping
            total_reads = self.reads_total
            raise StopIteration()
        else:
            return value


def calculate_mismatches(A, C, M):
    for name, vec_count in C.items():
        try:
            annots = A[name]
        except KeyError:
            continue
        counts = np.sum(np.sum(vec_count[0], axis=0) - np.amax(vec_count[0], axis=0)) + vec_count[1]
        M['Gene'][name] = [counts, vec_count[2]]
        for e, a in enumerate(annots):
            try:
                M[level_map[e]][a][0] += counts
                M[level_map[e]][a][1] += vec_count[2]
            except KeyError:
                M[level_map[e]][a] = [counts, vec_count[2]]


if __name__ == '__main__':
    D = {k: v for k, v in fasta_parse(sys.argv[1])}
    A = load_annotations(sys.argv[2])
    dataset = sys.argv[3]
    samplename = sys.argv[4]
    groupstring = sys.argv[5]
    for refname, start, cigar, reflen, readseq in SamParser('-'):
        # RefName, 0-start, CIGAR, RefLen, ReadSeq
        parse_cigar_reference(cigar, D[refname], start, refname, readseq)
    calculate_mismatches(A, cov_vectors, mismatches)
    with open(sys.argv[6], 'a') as out:
        for level, ldict in mismatches.items():
            for nodename, count in ldict.items():
                out.write('{},{},NA,NA,Alignment,{},{},{},{},{}\n'.format(dataset, samplename, groupstring, level,
                                                                          nodename, str(count[0]), str(count[1])))

