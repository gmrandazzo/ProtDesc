#!/usr/bin/env python3
import csv
def fasta2seqvect(ffasta):
    seq = []
    with open(ffasta, "r") as f:
        for line in f:
            if ".pdb" in line or ">" in line:
                continue
            else:
                seq.extend(line.strip())    
    return seq


def writedesc(l : list, fname, fout):
    row = []
    row = [fname.replace("fasta", "")]
    header = ["Protein"]
    for d in l:
        keys = list(d.keys())
        keys.sort()
        header.extend(keys)
        for key in keys:
            row.append(d[key])
    with open(fout, "w") as f:
        w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        w.writerow(header)
        w.writerow(row)
    return
