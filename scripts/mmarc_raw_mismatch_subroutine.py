#!/usr/bin/env python3


import sys
import os.path
import numpy as np


cov_vectors = {}
read_covs = {}
mismatches = {'Class': {}, 'Mechanism': {}, 'Group': {}, 'HMM': {}}
base_mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
level_map = {0: 'Class', 1: 'Mechanism', 2: 'Group'}

reads_mapped = 0
total_reads = 0
reads_not_considered = 0


def load_model_annotations(infile):
    ret = {}
    with open(infile, 'r') as f:
        data = f.read().split('\n')[1:]
        for line in data:
            if line:
                entry = line.split('\t')
                ret.setdefault(entry[1], (int(entry[3]), entry[4:]))
    return ret


class TbloutParser:
    def __init__(self, filepath):
        if os.path.exists(filepath):  # if file is a file, read from the file
            self.tblout_file = str(filepath)
            self.stdin = False
        elif not sys.stdin.isatty():  # else read from standard in
            self.stdin = True
        else:
            raise ValueError("Parameter filepath must be a tblout file")
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
                tblout_line = sys.stdin.readline()  # read from stdin
            else:
                tblout_line = self.tblout_file.readline()  # read from file
            if not tblout_line:
                return  # End of file
            if tblout_line[0] != "#":  # these lines are the actual reads
                self.reads_total += 1
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rReads processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = tblout_line.split()
                self.reads_mapping += 1
                return temp[0], temp[2], int(temp[4]), int(temp[5]), int(temp[6]), int(temp[7]), temp[11]
        self.tblout_file.close()  # catch all in case this line is reached
        assert False, "Should not reach this line"

    def __next__(self):
        if not self.stdin and type(self.tblout_file) is str:  # only open file here if sam_file is a str and not fileIO
            self.tblout_file = open(self.tblout_file, "r")
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.tblout_file.close()
            sys.stdout.write("\n{:d} reads mapped out of {:d} total reads\n".format(self.reads_mapping, self.reads_total))
            sys.stdout.flush()
            global reads_mapped
            global total_reads
            reads_mapped = self.reads_mapping
            total_reads = self.reads_total
            raise StopIteration()
        else:
            return value


def parse_tblout_reference(readseq, model_name, mstart, mstop, rstart, rstop, strand, model_len):
    global reads_mapped
    global reads_not_considered
    global read_covs
    global reads_seen
    if (mstop - mstart) != abs(rstop - rstart):
        reads_not_considered += 1
        return
    else:
        reads_mapped += 1
    if strand == '-':
        rstart, rstop = int(rstop), int(rstart)
    if model_name not in cov_vectors:
        cov_vectors[model_name] = [np.zeros((4, int(model_len)), dtype=np.float32), 0, 0]
    mstop -= 1
    mstart -= 1
    rstart -= 1
    rstop -= 1
    if readseq not in read_covs:
        read_covs[readseq] = ((model_name, mstart, mstop, rstart, rstop),)
    else:
        read_covs[readseq] += ((model_name, mstart, mstop, rstart, rstop),)
    cov_vectors[model_name][1] += 1
    for i in range(mstop - mstart):
        try:
            cov_vectors[model_name][0][base_mapping[readseq[rstart + i]]][mstart + i] += 1
        except KeyError:
            cov_vectors[model_name][2] += 1


def correct_multihits(C, R):
    mhits = {}
    for readseq, hits in R.items():
        if len(hits) > 1:
            for mname, mstart, mstop, rstart, rstop in hits:
                for i in range(mstop - mstart):
                    try:
                        C[mname][0][base_mapping[readseq[rstart + i]]][mstart + i] -= float(1) - (float(1) / len(hits))
                    except KeyError:
                        C[mname][2] -= float(1) - (float(1) / len(hits))
                mhits.setdefault(mname, []).append(float(1) - (float(1) / len(hits)))
    return mhits


def calculate_mismatches(A, C, M, H):
    for name, vec_count in C.items():
        annots = A[name][1]
        counts = np.sum(np.sum(vec_count[0], axis=0) - np.amax(vec_count[0], axis=0)) + vec_count[2]
        if name in H:
            vec_count[1] -= sum(H[name])
        for e, astring in enumerate(annots):
            annot_list = astring.split('|')
            for a in annot_list:
                try:
                    M[level_map[e]][a][0] += float(counts) / len(annot_list)
                    M[level_map[e]][a][1] += float(vec_count[1]) / len(annot_list)
                except KeyError:
                    M[level_map[e]][a] = [float(counts) / len(annot_list), float(vec_count[1]) / len(annot_list)]

if __name__ == '__main__':
    A = load_model_annotations(sys.argv[1])
    dataset = sys.argv[2]
    samplename = sys.argv[3]
    groupstring = sys.argv[4]
    for readseq, model_name, mstart, mstop, rstart, rstop, strand in TbloutParser('-'):
        # RefName, 0-start, CIGAR, RefLen, ReadSeq
        parse_tblout_reference(readseq, model_name, int(mstart), int(mstop), int(rstart), int(rstop), strand, A[model_name][0])
    multihits = correct_multihits(cov_vectors, read_covs)
    calculate_mismatches(A, cov_vectors, mismatches, multihits)
    sys.stdout.write('{}\t{}\t{}\n'.format(reads_not_considered, reads_mapped, total_reads))
    with open(sys.argv[5], 'a') as out:
        for level, ldict in mismatches.items():
            for nodename, count in ldict.items():
                out.write('{},{},NA,NA,MMARC,{},{},{},{},{}\n'.format(dataset, samplename, groupstring, level,
                                                                          nodename, str(count[0]), str(count[1])))

