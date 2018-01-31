#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict

from twobitreader import TwoBitFile

import pickle
import pickletools
import time

M = 7
Q = 13

# reference = TwoBitFile("GRCh38_no_alts.2bit")
reference = TwoBitFile("GRCh38_no_alts.2bit")

base_lookup = {
    'A': 0,
    'G': 1,
    'T': 2,
    'C': 3,
}


# Reference::::
# Bioinformatics. 2010 Sep 15; 26(18): i414–i419.  Published online 2010 Sep 4. doi:  10.1093/bioinformatics/btq364
# PMCID: PMC2935425
# A fast algorithm for exact sequence search in biological sequences using polyphase decomposition

def read_fasta(filename):
    with open(filename, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            print(record.id)
            return create_hash_table(record.seq._data)

def read_2bit(filename):
    return create_hash_table(get_dna_from2bit(filename))

def get_dna_from2bit(filename):
    return str(TwoBitFile(filename)['chr22'])

#        for record in SeqIO.parse(handle, 'fasta'):
#            print(record.id)
#            return create_hash_table(record.seq._data)


def create_hash_table(dna):
    table = defaultdict(list)
    reduced = down_sample(dna)
    k = None
    for i in range(len(reduced)-Q):
        key = reduced[i:i+Q]
        if "N" not in key:
            if k is None:
                k = encode(key)
            else:
                k += base_lookup[key[-1]] << Q
                k >>= 1

            table[k].append(i)
        else:
            k = None
    return table


def down_sample(dna, step=M):
    return dna[::step]


def encode(dna):
    # return dna
    val = 0
    for i, c in enumerate(dna.upper()):
        x = 0
        x = base_lookup[c]
        # if c == "A":
        #     x = 0
        # elif c == "C":
        #     x = 1
        # elif c == "G":
        #     x = 2
        # elif c == "T":
        #     x = 3
        # else:
        #     raise Exception("Unsupported char %s" % c)
        val += x << i
    return val


def match_dna(table, query):
    if len(query) < M*Q:
        raise Exception("Query must be at least %d characters" % M*Q)
    all_candidates = []
    for i in range(Q):
        key = encode(down_sample(query[i:])[:Q])
        if key in table:
            candidates = [(x - i) * M for x in table[key]]
            all_candidates += candidates

    print("Candidates:", candidates)
    return [c for c in all_candidates if check_candidate_match(c, query)]


def check_candidate_match(position, query):
    # seek to file
    ref = reference["chr22"][position:position+len(query)]
    #print("ref: {}\nqry: {}".format(ref, query))
    if ref == query:
        #print("They match")
        return True
    #print("No match")
    return False


def write_table_to(table, filename):
    with open(filename, 'wb') as handle:
        # pickletools.optimize(pickle.dumps(table))
        pickle.dump(table, handle)

def read_table_from(filename):
    with open(filename, 'rb') as handle:
        return pickle.load(handle)


def main():
    start_time = time.time()
    #table = read_fasta("chr1.fa")
    table = read_2bit("chr22.2bit")
    table_time = time.time()
    write_table_to(table, "chr22.index.pickle")
    write_table_time = time.time()
    table = read_table_from("chr22.index.pickle")
    read_table_time = time.time()
    # query = "CCACCTGTACATGCTATCTGAAGGACAGCCTCCAGGGCACACAGAGGATGGTATTTACACATGCACACATGGCTACTGATGGGGCAAGCACTTCACAACCCCTCATGATCACGTGCAGCAGACAATGTGGCCTCTGCAGAGGGGGAACGGAGACCGGAGGCTGAGACTGGCAAGGCTGGACCTGAGTGTCGTCACCTAAATTCAGACGGG"
    query = "GTAATCTTAGCACTTTGGGAGGCGGAGACGGATGTATCGCTTGAGCTCAGGAGTTGAAGACCAGCCTGGGCAACATACTGAGACTCCGTCTTGTATAATTTAATTAAAATTTAAAAAAAGAAGAGAAAAAGACCTGTGTT"

    matchCount = 1
    matches = []
    for i in range(matchCount):
        matches = match_dna(table, query)
        # print(matches)
    match_time = time.time()
    print(matches)

    print("\n".join(
          ["Elapsed time:",
           "Table creation: {}".format(table_time - start_time),
           "Table writing:  {}".format(write_table_time - table_time),
           "Table reading:  {}".format(read_table_time - write_table_time),
           "Per query:    : {}".format((match_time - read_table_time) / matchCount),
           ]))


main()
