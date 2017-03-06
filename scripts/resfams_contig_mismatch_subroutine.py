#!/usr/bin/env python3


import sys
import os.path
import re
import numpy as np


cov_vectors = {}
mismatches = {}
base_mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
level_map = {0: 'Class', 1: 'Mechanism', 2: 'Group'}

reads_mapped = 0
total_reads = 0


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
            ali = (int(entry[19]), int(entry[20]))
            locs[0] += 3 * (ali[0] - 1)
            locs[0] += 3 * (ali[1] - 1)
            ret.setdefault(contig, []).append((int(locs[0]), int(locs[1]), entry[4]))
    return ret


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
            if sam_line[:3] == '@SQ':
                entry = sam_line.split('\t')[1:]
                name = entry[0].replace('SN:', '')
                length = int(entry[1].replace('LN:', ''))
                self.header_lens[name] = length
            elif sam_line[0] != '@':  # these lines are the actual reads
                self.reads_total += 1
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rReads processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = sam_line.split()
                if (int(temp[1]) & 4) == 0:
                    self.reads_mapping += 1
                    return temp[2], int(temp[3])-1, temp[5], self.header_lens[temp[2]], temp[9]  # RefName, 0-start, CIGAR, RefLen, ReadSeq
            else:
                continue
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


def parse_cigar_reference(s, ref_len, ref_start, ref_name, readseq):
    ret = re.findall(r'(\d+)([A-Z=]{1})', s)
    if ref_name not in cov_vectors:
        cov_vectors[ref_name] = [np.zeros((4, ref_len), dtype=np.int32), 0, 0, tuple()]
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
                #print(ref_start, ref_len, current_loc, occ, op, s)
                try:
                    cov_vectors[ref_name][0][base_mapping[readseq[i]]][current_loc] += 1
                    current_loc += 1
                except (KeyError, IndexError) as e:
                    if i >= len(readseq) or current_loc >= ref_len:
                        print(ref_start, ref_len, current_loc, occ, op, s, readseq)
                    current_loc += 1
        elif op in {'S', 'H'}:
            pass
    cov_vectors[ref_name][3] += ((ref_start, current_loc),)


def calculate_mismatches(A, D, C, M):
    for name, tdict in D.items():
        try:
            lvec = C[name][0]
        except KeyError:
            continue
        readdict = {}
        contigcov = np.zeros((1, lvec.shape[1]))
        for triplets in tdict:
            for i in range(triplets[1] - triplets[0]):
                contigcov[0, triplets[0] + i] += 1
            for e, refpair in enumerate(C[name][3]):
                if max(refpair[0], triplets[0]) <= min(refpair[1], triplets[1]):
                    readdict.setdefault(e, []).append(triplets[2])
        contigcov[contigcov == 0] = 1
        countvec = (np.sum(lvec, axis=0) - np.amax(lvec, axis=0)) / contigcov
        for triplets in tdict:
            annot = A[triplets[2]]
            counts = np.sum(countvec[0, triplets[0]:triplets[1]])
            try:
                M[annot][0] += float(counts)
            except KeyError:
                M[annot] = [float(counts), 0]
        for readnum, modelnames in readdict.items():
            for mname in modelnames:
                annot = A[mname]
                M[annot][1] += float(1) / len(modelnames)


if __name__ == '__main__':
    A = load_resfams_metadata(sys.argv[1])
    dataset = sys.argv[2]
    samplename = sys.argv[3]
    groupstring = sys.argv[4]
    for refname, start, cigar, reflen, readseq in SamParser('-'):
        # RefName, 0-start, CIGAR, RefLen, ReadSeq
        parse_cigar_reference(cigar, reflen, start, refname, readseq)
    D = load_domains(sys.argv[5])
    calculate_mismatches(A, D, cov_vectors, mismatches)
    with open(sys.argv[6], 'a') as out:
        for nodename, count in mismatches.items():
            out.write('{},{},NA,NA,Resfams,{},Class,{},{},{}\n'.format(dataset, samplename, groupstring,
                                                                          nodename, str(count[0]), str(count[1])))

