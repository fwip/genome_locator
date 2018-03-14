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

import gzip
import pickle
import pickletools
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

    def read_2bit(self, filename):
        table = {}
        tbf = TwoBitFile(filename)
        for (chrom, dna) in tbf.items():
            print("Chrom:", chrom, self.elapsed())
            table[chrom] = self.get_table_for_chrom(filename, chrom)

        return table

    def get_table_for_chrom(self, filename, chrom, offset):
        tbf = TwoBitFile(filename)
        table = self.create_hash_table(tbf[chrom], offset)
        lh = LookupHash(table)
        print("Converted to lookuphash", self.elapsed())
        return lh

    def write_tables_from_2bit(self, filename, outname=None):
        # tbf = TwoBitFile(filename)
        if outname is None:
            outname = "{}.M{}.Q{}.index.hdf5".format(filename, self.M, self.Q)
        # for (chrom, dna) in tbf.items():
        offset = 0
        master_table = defaultdict(list)
        for (chrom, length) in self.chroms:
            self.create_hash_table(TwoBitFile(filename)[chrom], table=master_table, offset=offset)
            print("Updated master table", self.elapsed())
            offset += length
            # outfile = "{}.M{}.Q{}.index.gz".format(chrom, M, Q)
            # write_table_to(table, outfile)
            # write_chrom_h5(table, chrom, h5_file)
        lh = LookupHash(table=master_table)
        print("Created lookup hash", self.elapsed())
        h5_file = h5py.File(outname, "w")
        self.write_h5(lh, h5_file)
        h5_file.close()
        print("Wrote h5 file", self.elapsed())

    def rotate_key(self, k, new_letter):
        k += (self.base_lookup[new_letter] << (self.Q * 2))
        return k >> 2

    def create_hash_table(self, dna, table=None, offset=0):
        if table is None:
            table = defaultdict(list)

        dna = str(dna)
        print("Read dna", self.elapsed())
        reduced = self.down_sample(dna)
        k = None
        for i in range(len(reduced) - self.Q + 1):
            key = reduced[i:i + self.Q]
            if "N" not in key:
                if k is None:
                    k = self.encode(key)
                else:
                    k = self.rotate_key(k, key[-1])

                table[k].append(i+offset)
            else:
                k = None
        print("Created table", self.elapsed())
        return table

    def down_sample(self, dna, step=M):
        return dna[::step]

    def encode(self, dna):
        val = 0
        for i, c in enumerate(dna.upper()):
            x = self.base_lookup[c]
            val += x << (i * 2)
        return val

    def match_file(self, filename, query):
        f = h5py.File(filename)
        #chroms = f['chromosomes']
        table = LookupHash(h5_group=f['index'])
        #table = {
        #    chrom: LookupHash(h5_group=chroms[chrom])
        #    for chrom in chroms.keys()
        #}
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
            key = self.encode(self.down_sample(query[i:])[:self.Q])
            if key in table:
                for x in table[key]:
                    chrom, deindexed = self.deindex(x)
                    all_candidates.append((chrom, self.M*deindexed - i))

        print("Candidates", all_candidates)
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
        print("r==q?", ref, query)
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

    def write_chrom_h5(self, table, chrom_name, h5file):
        group = h5file.require_group("chromosomes").create_group(chrom_name)
        group.create_dataset("offsets", compression="gzip",
                             shuffle=True, data=table.offsets)
        group.create_dataset("counts", compression="gzip",
                             shuffle=True, data=table.counts)
        group.create_dataset("positions", compression="gzip",
                             shuffle=True, data=table.positions)
        print("Wrote hdf5", self.elapsed())

    def write_table_to(self, table, filename):
        print("Writing to", filename)
        with gzip.open(filename, 'wb') as handle:
            dump = pickle.dumps(table)
            print("Dumped", self.elapsed())
            optimized = pickletools.optimize(dump)
            print("Optimized", self.elapsed())
            handle.write(optimized)
            print("Wrote", self.elapsed())

    def read_table_from_h5(self, filename):
        chroms = h5py.File(filename)['chromosomes']
        return {
            chrom: LookupHash(h5_group=chroms[chrom])
            for chrom in chroms.keys()
        }

    def read_table_from(self, filename):
        with gzip.open(filename, 'rb') as handle:
            return pickle.load(handle)



def main(*args):

    search = GenomeSearcher()
    search.M = 10
    search.Q = 10

    read_table_time = time.time()
    query = "GTAATCTTAGCACTTTGGGAGGCGGAGACGGATGTATCGCTTGAGCTCAGGAGTTGAAGACCAGCCTGGGCAACATACTGAGACTCCGTCTTGTATAATTTAATTAAAATTTAAAAAAAGAAGAGAAAAAGACCTGTGTT"
    index_file = "GRCh38_no_alts.2bit.M{}.Q{}.index.hdf5".format(search.M, search.Q)
    if len(args) > 1:
        print("Args", args)
        query = args[1]
    if len(args) > 2:
        index_file = args[2]

    matchCount = 50
    matches = []
    for i in range(matchCount):
        # matches = match_dna(table, query)
        matches = search.match_file(index_file, query)
    match_time = time.time()
    print(matches)

    print("\n".join(
          ["Elapsed time:",
           # "Table creation: {}".format(table_time - start_time),
           # "Table writing:  {}".format(write_table_time - table_time),
           # "Table reading:  {}".format(read_table_time - write_table_time),
           "Per query:    : {}".format(
               (match_time - read_table_time) / matchCount),
           ]))


if __name__ == "__main__":
    import sys
    main(sys.argv)
