#!/usr/bin/env python3

import sys
import os.path
import re

read_mapped = 0
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
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rReads processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = sam_line.split()
                if (int(temp[1]) & 4) == 0:
                    self.reads_mapping += 1
                    return temp[2], int(temp[3]), temp[5]  # RefName, 1-start, CIGAR
        self.sam_file.close()  # catch all in case this line is reached
        assert False, "Should not reach this line"

    def __next__(self):
        if not self.stdin and type(self.sam_file) is str:  # only open file here if sam_file is a str and not fileIO
            self.sam_file = open(self.sam_file, "r")
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.sam_file.close()
            sys.stdout.write(
                "\n{:d} reads mapped out of {:d} total reads\n".format(self.reads_mapping, self.reads_total))
            sys.stdout.flush()
            global reads_mapped
            global total_reads
            reads_mapped = self.reads_mapping
            total_reads = self.reads_total
            raise StopIteration()
        else:
            return value


def parse_sam(infile, D, R, dantas_string, dataset, run_mode):
    counts = {}
    with open(infile, 'r') as f:
        data = f.read().split('\n')
        for line in data:
            if line.startswith('@') or not line:
                continue
            entry = line.split('\t')
            model_hits = set()
            if int(entry[1]) & 4 == 0:
                seqlen, mismatch = parse_cigar(entry[5])
                contig = entry[2]
                locs = (int(entry[3]), int(entry[3]) + seqlen - 1)
                try:
                    for triplets in D[contig]:
                        if max(locs[0], triplets[0]) <= min(locs[1], triplets[1]):
                            if R[triplets[2]] != 'NA':
                                model_hits.add(triplets[2])
                except KeyError:
                    continue
            if model_hits:
                for x in model_hits:
                    if run_mode == 'mismatch' and mismatch != 0:
                        try:
                            counts[R[x]] += mismatch
                        except KeyError:
                            counts[R[x]] = mismatch
                    else:
                        try:
                            counts[R[x]] += float(1) / len(model_hits)
                        except KeyError:
                            counts[R[x]] = float(1) / len(model_hits)
    for k, v in counts.items():
        sys.stdout.write('{},{},Resfams,NA,Class,{},{}\n'.format(dataset, dantas_string, k, str(v)))


if __name__ == '__main__':
    R = load_resfams_metadata(sys.argv[1])
    D = load_domains(sys.argv[2])
    dantas_string = sys.argv[3]
    dataset = sys.argv[4]

    counts = {}
    for refname, start, cigar in SamParser('-'):
        # RefName, 1-start, CIGAR, RefLen, ReadSeq
        stop = start + parse_cigar(cigar) - 1
        if refname not in D:
            continue
        model_hits = set()
        for triplets in D[refname]:
            if max(start, triplets[0]) <= min(stop, triplets[1]):
                if R[triplets[2]] != 'NA':
                    model_hits.add(triplets[2])
        if model_hits:
            for x in model_hits:
                try:
                    counts[R[x]] += float(1) / len(model_hits)
                except KeyError:
                    counts[R[x]] = float(1) / len(model_hits)
    with open(sys.argv[5], 'a') as out:
        for k, v in counts.items():
            out.write('{},{},Resfams,NA,Class,{},{}\n'.format(dataset, dantas_string, k, str(v)))




