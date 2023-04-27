#!/usr/bin/env python3
import sys
from protinp import fasta2seqvect
from protinp import writedesc
from aminoengine import calc_mw
from aminoengine import triad_analyze
from aminoengine import quadriad_analyze
from aminoengine import iad_db_load
from aminoengine import iad_match

def main():
    if len(sys.argv) != 4:
        print(f"\nUsage {sys.argv[0]} [FASTA] [DB] [DESC OUT]")
        exit()

    seq = fasta2seqvect(sys.argv[1])
    # mw = calc_mw(seq)
    iads = triad_analyze(seq)
    # iads.extend(quadriad_analyze(seq))
    db = iad_db_load(sys.argv[2])
    d = iad_match(iads, db)
    writedesc([d], sys.argv[1], sys.argv[3])
    return

if __name__ in "__main__":
    main()

