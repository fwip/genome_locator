#!/usr/bin/env python3

import sys
import search


def main(infile, M, Q):
    search.M = M
    search.Q = Q

    search.write_tables_from_2bit(infile)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: build.py infile M Q")
        sys.exit(2)
    infile = sys.argv[1]
    M = int(sys.argv[2])
    Q = int(sys.argv[3])
    main(infile, M, Q)
