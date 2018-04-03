#!/usr/bin/env python3

# Reference::::
# Bioinformatics. 2010 Sep 15; 26(18): i414-i419.
# Published online 2010 Sep 4.
# doi:  10.1093/bioinformatics/btq364
# PMCID: PMC2935425
# A fast algorithm for exact sequence search in biological sequences using
# polyphase decomposition


from collections import defaultdict
from twobitreader import TwoBitFile
import h5py
import numpy
import time

from lookup_hash import LookupHash


class GenomeSearcher():

    M = 10
    Q = 10

    genome_name = "GRCh38_no_alts"
    reference = TwoBitFile("%s.2bit" % genome_name)

    chroms = []

    base_lookup = {
        'A': 0,
        'G': 1,
        'T': 2,
        'C': 3,
    }

    last_time = time.time()

    def __init__(self, **kwargs):
        if 'reference_file' in kwargs:
            self.reference = TwoBitFile(kwargs['reference_file'])

        if 'M' in kwargs:
            self.M = kwargs['M']
        if 'Q' in kwargs:
            self.Q = kwargs['Q']

        self.chroms = [
            (c, ((x + self.M-1)//self.M) - self.Q + 1)
            for (c, x) in sorted(self.reference.sequence_sizes().items(),
                                 key=lambda x: x[1],
                                 reverse=True)
        ]

    def elapsed(self):
        new_time = time.time()
        diff = new_time - self.last_time
        self.last_time = new_time
        return diff

    def write_tables_from_2bit(self, filename, outname=None):
        if outname is None:
            outname = "{}.M{}.Q{}.index.hdf5".format(filename, self.M, self.Q)
        offset = 0
        master_table = defaultdict(list)
        master_table = [None for _ in range(4**self.Q)]
        master_table = numpy.zeros(
            (sum(l for _, l in self.chroms), 2),
            dtype=numpy.uint64)
        tbl_len = 0
        for (chrom, length) in self.chroms:
            sub_table = master_table[tbl_len:]
            tbl_len += self.create_hash_table(
                TwoBitFile(filename)[chrom],
                table=sub_table,
                offset=offset,
                length=length)
            print("Updated master table", self.elapsed())
            offset += length
        lh = LookupHash(table=master_table)
        print("Created lookup hash", self.elapsed())
        h5_file = h5py.File(outname, "w")
        self.write_h5(lh, h5_file)
        h5_file.close()
        print("Wrote h5 file", self.elapsed())

    def rotate_key(self, k, new_letter):
        k += (self.base_lookup[new_letter] << (self.Q * 2))
        return k >> 2

    def create_hash_table(self, dna, table=None, offset=0, length=0):
        if table is None:
            table = defaultdict(list)

        dna = str(dna)
        print("Read dna", self.elapsed())
        reduced = self.down_sample(dna)
        values = 0
        k = None
        for i in range(len(reduced) - self.Q + 1):
            key = reduced[i:i + self.Q]
            if "N" not in key:
                if k is None:
                    k = self.encode(key)
                else:
                    k = self.rotate_key(k, key[-1])

                table[values] = (k, (i+offset))
                values += 1
            else:
                k = None
        print("Created table", self.elapsed())
        return values

    def down_sample(self, dna):
        return dna[::self.M]

    def encode(self, dna):
        val = 0
        for i, c in enumerate(dna.upper()):
            x = self.base_lookup[c]
            val += x << (i * 2)
        return val

    def match_file(self, filename, query):
        f = h5py.File(filename)
        table = LookupHash(h5_group=f['index'])
        results = self.match_dna(table, query)
        f.close()
        return results

    def deindex(self, idx):
        for i in range(0, len(self.chroms)):
            c, length = self.chroms[i]
            if idx < length:
                return (c, idx)
            idx -= length
        raise Exception("Something is wrong")

    def match_dna(self, table, query):
        if len(query) < self.M * self.Q:
            raise Exception("Query must be at least %d characters" % (self.M * self.Q))
        all_candidates = []
        for i in range(self.M):
            long_sample = query[i:]
            down_sample = self.down_sample(long_sample)
            substr = self.down_sample(query[i:])[:self.Q]
            key = self.encode(substr)
            if key in table:
                for x in table[key]:
                    chrom, deindexed = self.deindex(x)
                    all_candidates.append((chrom, self.M*deindexed - i))

        return [
            (c[0], c[1] + 1)
            for c in all_candidates
            if self.check_candidate_match(c, query)
        ]

    def check_candidate_match(self, position, query):
        # seek to file
        chrom = position[0]
        idx = position[1]
        ref = self.reference[chrom][idx:idx + len(query)]
        return ref == query

    def write_h5(self, table, h5file):
        group = h5file.require_group("index")
        group.create_dataset("offsets", compression="gzip",
                             shuffle=True, data=table.offsets)
        group.create_dataset("counts", compression="gzip",
                             shuffle=True, data=table.counts)
        group.create_dataset("positions", compression="gzip",
                             shuffle=True, data=table.positions)
        print("Wrote hdf5", self.elapsed())


def main(*args):
    M = 10
    Q = 10
    index_file = "GRCh38_no_alts.2bit.M{}.Q{}.index.hdf5".format(M, Q)
    if len(args) > 1:
        print("Args", args)
        query = args[1]
    if len(args) > 2:
        index_file = args[2]
    if len(args) > 3:
        M = int(args[3])
    if len(args) > 4:
        Q = int(args[4])
    search = GenomeSearcher(M=M, Q=Q)

    query = "GATTTGACTTTACCTTGAGCTTTGTCAGTTTACGATGCTATTTCAGTTTTGTGCTCAGATTTGAGTGATTGCAGGAAGAGAATAAATTTCTTTAATGCTGTCAAGACTTTAAATAGATACAGACAGAGCATTTTCACTTTTTCCTACATC"

    matchCount = 10
    matches = []
    print("Matches:", search.match_file(index_file, query))
    init_time = time.time()
    for i in range(matchCount):
        search.match_file(index_file, query)
    match_time = time.time()
    print(matches)

    print("\n".join(
          ["Elapsed time:",
           "Per query:    : {}".format(
               (match_time - init_time) / matchCount),
           ]))


if __name__ == "__main__":
    import sys
    main(sys.argv)
