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


def elapsed():
    global last_time
    new_time = time.time()
    diff = new_time - last_time
    last_time = new_time
    return diff


def read_2bit(filename):
    table = {}
    tbf = TwoBitFile(filename)
    for (chrom, dna) in tbf.items():
        print("Chrom:", chrom, elapsed())
        table[chrom] = get_table_for_chrom(filename, chrom)

    return table


def get_table_for_chrom(filename, chrom, offset):
    tbf = TwoBitFile(filename)
    table = create_hash_table(tbf[chrom], offset)
    lh = LookupHash(table)
    print("Converted to lookuphash", elapsed())
    return lh


def write_tables_from_2bit(filename, outname=None):
    init()
    # tbf = TwoBitFile(filename)
    if outname is None:
        outname = "{}.M{}.Q{}.index.hdf5".format(filename, M, Q)
    # for (chrom, dna) in tbf.items():
    offset = 0
    master_table = defaultdict(list)
    for (chrom, length) in chroms:
        create_hash_table(TwoBitFile(filename)[chrom], table=master_table, offset=offset)
        print("Updated master table", elapsed())
        offset += length
        # outfile = "{}.M{}.Q{}.index.gz".format(chrom, M, Q)
        # write_table_to(table, outfile)
        # write_chrom_h5(table, chrom, h5_file)
    lh = LookupHash(table=master_table)
    print("Created lookup hash", elapsed())
    h5_file = h5py.File(outname, "w")
    write_h5(lh, h5_file)
    h5_file.close()
    print("Wrote h5 file", elapsed())


def rotate_key(k, new_letter):
    k += (base_lookup[new_letter] << (Q * 2))
    return k >> 2


def create_hash_table(dna, table=None, offset=0):
    if table is None:
        table = defaultdict(list)

    dna = str(dna)
    print("Read dna", elapsed())
    reduced = down_sample(dna)
    k = None
    for i in range(len(reduced) - Q + 1):
        key = reduced[i:i + Q]
        if "N" not in key:
            if k is None:
                k = encode(key)
            else:
                k = rotate_key(k, key[-1])

            table[k].append(i+offset)
        else:
            k = None
    print("Created table", elapsed())
    return table


def down_sample(dna, step=M):
    return dna[::step]


def encode(dna):
    val = 0
    for i, c in enumerate(dna.upper()):
        x = base_lookup[c]
        val += x << (i * 2)
    return val


def match_file(filename, query):
    init()
    f = h5py.File(filename)
    #chroms = f['chromosomes']
    table = LookupHash(h5_group=f['index'])
    #table = {
    #    chrom: LookupHash(h5_group=chroms[chrom])
    #    for chrom in chroms.keys()
    #}
    results = match_dna(table, query)
    f.close()
    return results


def deindex(idx):
    for i in range(0, len(chroms)):
        c, length = chroms[i]
        if idx < length:
            return (c, idx)
        idx -= length
    raise Exception("Something is wrong")


def match_dna(table, query):
    if len(query) < M * Q:
        raise Exception("Query must be at least %d characters" % M * Q)
    all_candidates = []
    for i in range(M):
        key = encode(down_sample(query[i:])[:Q])
        if key in table:
            for x in table[key]:
                chrom, deindexed = deindex(x)
                all_candidates.append((chrom, M*deindexed - i))

    return [
        (c[0], c[1] + 1)
        for c in all_candidates
        if check_candidate_match(c, query)
    ]


def check_candidate_match(position, query):
    # seek to file
    chrom = position[0]
    idx = position[1]
    ref = reference[chrom][idx:idx + len(query)]
    return ref == query


def write_h5(table, h5file):
    group = h5file.require_group("index")
    group.create_dataset("offsets", compression="gzip",
                         shuffle=True, data=table.offsets)
    group.create_dataset("counts", compression="gzip",
                         shuffle=True, data=table.counts)
    group.create_dataset("positions", compression="gzip",
                         shuffle=True, data=table.positions)
    print("Wrote hdf5", elapsed())


def write_chrom_h5(table, chrom_name, h5file):
    group = h5file.require_group("chromosomes").create_group(chrom_name)
    group.create_dataset("offsets", compression="gzip",
                         shuffle=True, data=table.offsets)
    group.create_dataset("counts", compression="gzip",
                         shuffle=True, data=table.counts)
    group.create_dataset("positions", compression="gzip",
                         shuffle=True, data=table.positions)
    print("Wrote hdf5", elapsed())


def write_table_to(table, filename):
    print("Writing to", filename)
    with gzip.open(filename, 'wb') as handle:
        dump = pickle.dumps(table)
        print("Dumped", elapsed())
        optimized = pickletools.optimize(dump)
        print("Optimized", elapsed())
        handle.write(optimized)
        print("Wrote", elapsed())


def read_table_from_h5(filename):
    chroms = h5py.File(filename)['chromosomes']
    return {
        chrom: LookupHash(h5_group=chroms[chrom])
        for chrom in chroms.keys()
    }


def read_table_from(filename):
    with gzip.open(filename, 'rb') as handle:
        return pickle.load(handle)


def init():
    global chroms
    chroms = [
        (c, ((x + M-1)//M) - Q + 1)
        for (c, x) in sorted(reference.sequence_sizes().items(),
                             key=lambda x: x[1],
                             reverse=True)
    ]
    print("Chroms", chroms)


def main():
    init()
    start_time = time.time()
    #table = read_2bit("%s.2bit" % genome_name)
    table_time = time.time()
    #write_table_to(table, "%s.index.pickle" % genome_name)
    write_table_time = time.time()
    #table = read_table_from("%s.2bit.M10.Q10.index.hdf5" % genome_name)
    read_table_time = time.time()
    query = "GTAATCTTAGCACTTTGGGAGGCGGAGACGGATGTATCGCTTGAGCTCAGGAGTTGAAGACCAGCCTGGGCAACATACTGAGACTCCGTCTTGTATAATTTAATTAAAATTTAAAAAAAGAAGAGAAAAAGACCTGTGTT"

    matchCount = 50
    matches = []
    global M
    global Q
    M = 10
    Q = 10
    for i in range(matchCount):
        # matches = match_dna(table, query)
        matches = match_file("GRCh38_no_alts.2bit.M{}.Q{}.index.hdf5".format(M, Q), query)
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
    main()
