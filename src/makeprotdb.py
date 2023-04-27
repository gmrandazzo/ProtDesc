#!/usr/bin/env python3
import sys
from protinp import fasta2seqvect
from aminoengine import triad_analyze
from aminoengine import quadriad_analyze
from aminoengine import iad_db_save

def main():
    if len(sys.argv) != 3:
        print(f"\nUsage {sys.argv[0]} [FASTA] [DB out]")
        exit()
        
    seq = fasta2seqvect(sys.argv[1])
    iads = triad_analyze(seq)
    # iads.extend(quadriad_analyze(seq))
    iad_db_save(iads, sys.argv[2])
    return

if __name__ in "__main__":
    main()

