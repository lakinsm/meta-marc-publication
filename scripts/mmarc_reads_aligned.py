#!/usr/bin/env python3

import re
import os.path
import sys

ltrans = {1: 'Class', 2: 'Mechanism', 3: 'Group'}
counts = {1: {}, 2: {}, 3: {}}

reads_mapped = 0
total_reads = 0


def load_mmarc_metadata(infile):
    ret = {}
    with open(infile, 'r') as f:
        data = f.read().split('\n')[1:]
        for line in data:
            if not line:
                continue
            entry = line.split('\t')
            ret[entry[1]] = entry[-3:]
    return ret


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


if __name__ == '__main__':
    R = load_mmarc_metadata(sys.argv[1])
    D = load_tblout(sys.argv[2])
    dantas_string = sys.argv[3]
    dataset = sys.argv[4]
    groupstring = sys.argv[5]
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
                for e, a in enumerate(R[x]):
                    annots = a.split('|')
                    correct = len(model_hits) * len(annots)
                    for annot in annots:
                        try:
                            counts[e + 1][annot] += float(1) / correct
                        except KeyError:
                            counts[e + 1][annot] = float(1) / correct
    with open(sys.argv[6], 'a') as out:
        for level, ldict in counts.items():
            for k, count in ldict.items():
                out.write(
                    '{},{},MMARC,{},{},{},{}\n'.format(dataset, dantas_string, groupstring, ltrans[level], k,
                                                       str(count)))

