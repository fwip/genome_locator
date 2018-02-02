#!/usr/bin/env python3

import sys
# from search import read_2bit, write_table_to, M, Q
import search


def main(infile, M, Q, outputfile):
    search.M = M
    search.Q = Q

    table = search.read_2bit(infile)
    search.write_table_to(table, outputfile)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: build.py infile M Q")
    infile = sys.argv[1]
    M = int(sys.argv[2])
    Q = int(sys.argv[3])
    outfile = "{root}.M{M}.Q{Q}.index".format(root=infile, M=M, Q=Q)
    main(infile, M, Q, outfile)
