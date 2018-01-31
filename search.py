#!/usr/bin/env python3

# Reference::::
# Bioinformatics. 2010 Sep 15; 26(18): i414–i419.  Published online 2010 Sep 4. doi:  10.1093/bioinformatics/btq364
# PMCID: PMC2935425
# A fast algorithm for exact sequence search in biological sequences using
# polyphase decomposition


from collections import defaultdict
from twobitreader import TwoBitFile

import gzip
import pickle
import pickletools
import time

M = 7
Q = 13

genome_name = "GRCh38_no_alts"
reference = TwoBitFile("%s.2bit" % genome_name)

base_lookup = {
    'A': 0,
    'G': 1,
    'T': 2,
    'C': 3,
}


def read_fasta(filename):
    with open(filename, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            print(record.id)
            return create_hash_table(record.seq._data)


def read_2bit(filename):
    tbf = TwoBitFile(filename)
    table = {}
    for (chrom, dna) in tbf.items():
        chrom_table = create_hash_table(str(dna))
        table[chrom] = chrom_table

    return table


def rotate_key(k, new_letter):
    k += (base_lookup[new_letter] << (Q * 2))
    return k >> 2


def create_hash_table(dna):
    table = defaultdict(list)
    reduced = down_sample(dna)
    k = None
    for i in range(len(reduced) - Q + 1):
        key = reduced[i:i + Q]
        if "N" not in key:
            if k is None:
                k = encode(key)
            else:
                k = rotate_key(k, key[-1])

            table[k].append(i)
        else:
            k = None
    return table


def down_sample(dna, step=M):
    return dna[::step]


def encode(dna):
    val = 0
    for i, c in enumerate(dna.upper()):
        x = base_lookup[c]
        val += x << (i * 2)
    return val


def match_dna(table, query):
    if len(query) < M * Q:
        raise Exception("Query must be at least %d characters" % M * Q)
    all_candidates = []
    for i in range(Q):
        key = encode(down_sample(query[i:])[:Q])
        for (chrom, subtable) in table.items():
            if key in subtable:
                candidates = [(chrom, (x - i) * M) for x in subtable[key]]
                all_candidates += candidates

    return [c for c in all_candidates if check_candidate_match(c, query)]


def check_candidate_match(position, query):
    # seek to file
    chrom = position[0]
    idx = position[1]
    ref = reference[chrom][idx:idx + len(query)]
    return ref == query


def write_table_to(table, filename):
    with gzip.open(filename, 'wb') as handle:
        # handle.write(pickletools.optimize(pickle.dumps(table)))
        pickle.dump(table, handle)


def read_table_from(filename):
    with gzip.open(filename, 'rb') as handle:
        return pickle.load(handle)


def main():
    start_time = time.time()
    table = read_2bit("%s.2bit" % genome_name)
    table_time = time.time()
    write_table_to(table, "%s.index.pickle.gz" % genome_name)
    write_table_time = time.time()
    table = read_table_from("%s.index.pickle.gz" % genome_name)
    read_table_time = time.time()
    query = "GTAATCTTAGCACTTTGGGAGGCGGAGACGGATGTATCGCTTGAGCTCAGGAGTTGAAGACCAGCCTGGGCAACATACTGAGACTCCGTCTTGTATAATTTAATTAAAATTTAAAAAAAGAAGAGAAAAAGACCTGTGTT"

    matchCount = 1
    matches = []
    for i in range(matchCount):
        matches = match_dna(table, query)
        print(matches)
    match_time = time.time()
    print(matches)

    print("\n".join(
          ["Elapsed time:",
           "Table creation: {}".format(table_time - start_time),
           "Table writing:  {}".format(write_table_time - table_time),
           "Table reading:  {}".format(read_table_time - write_table_time),
           "Per query:    : {}".format(
               (match_time - read_table_time) / matchCount),
           ]))


main()
