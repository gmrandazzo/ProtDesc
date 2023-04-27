#!/usr/bin/env python3
#
# List of aminoacids
#
# Ala 	A 	Alanine
# Arg 	R 	Arginine
# Asn 	N 	Asparagine
# Asp 	D 	Aspartic acid
# Cys 	C 	Cysteine
# Gln 	Q 	Glutamine
# Glu 	E 	Glutamic acid
# Gly 	G 	Glycine
# His 	H 	Histidine
# Ile 	I 	Isoleucine
# Leu 	L 	Leucine
# Lys 	K 	Lysine
# Met 	M 	Methionine
# Phe 	F 	Phenylalanine
# Pro 	P 	Proline
# Pyl 	O 	Pyrrolysine
# Ser 	S 	Serine
# Sec 	U 	Selenocysteine
# Thr 	T 	Threonine
# Trp 	W 	Tryptophan
# Tyr 	Y 	Tyrosine
# Val 	V 	Valine
#
# Asx 	B 	Aspartic acid or Asparagine
# Glx 	Z 	Glutamic acid or Glutamine
# Xaa 	X 	Any amino acid
# Xle 	J 	Leucine or Isoleucine
#

from copy import copy

def general_iad_analyzer(seq : list, length : int):
    iads = []
    for i in range(0, len(seq), 1):
        iad = seq[i:i+length]
        if len(iad) == length:
            iads.append(iad)
        else:
            break
    return iads

def quadriad_analyze(seq : list):
    return general_iad_analyzer(seq, 4)

def triad_analyze(seq : list):
    return general_iad_analyzer(seq, 3)

def iad2str(triad):
    return "".join(triad)

def iad_db_save(triads, fout):
    predb = iad_db_load(fout)
    lst = []
    for triad in triads:
        str_triad = iad2str(triad)
        if str_triad in predb.keys():
            continue    
        else:
            lst.append(str_triad)
    lst = list(set(lst))
    with open(fout, "a") as f:
        for item in lst:
            f.write(f'{item}\n')
    return 0

def iad_db_load(fdb):
    db = {}
    try:
        with open(fdb, "r") as f:
            for line in f:
                db[line.strip()] = 0
    except:
        print("DB not found")
    return db

def iad_match(iads, db):
    copy_db = copy(db)
    for iad in iads:
        str_iad = iad2str(iad)
        try:
            copy_db[str_iad] += 1
        except:
            continue
    return copy_db


def calc_mw(seq : list):
    """
    Calculate the mol weight of the protein
    and the single aminoacid molweight sum
    """
    d = {}
    d['A'] = 89.094
    d['C'] = 121.15
    d['D'] = 133.103
    d['E'] = 147.130
    d['F'] = 165.192
    d['G'] = 75.067
    d['H'] = 155.157 
    d['I'] = 131.175
    d['K'] = 146.190
    d['L'] = 131.175
    d['M'] = 149.21
    d['N'] = 132.119
    d['O'] = 255.313 
    d['P'] =  115.132
    d['Q'] = 146.146
    d['R'] = 174.204
    d['S'] = 105.093
    d['T'] = 119.120 
    d['U'] = 168.065
    d['V'] = 117.148
    d['W'] = 204.229
    d['Y'] = 181.191
    r = {}
    for key in d.keys():
        r[f'MW_{key}'] = 0.0
    n = 0
    mw = 0.0
    for ac in seq:
        try:
            r[f'MW_{ac}'] += d[ac]
            n += 1
            mw += d[ac]
        except:
            continue
    
    # Subtract water because 
    # on every peptide bond 1 H2O is
    # created and expelled
    mw -= (18.0153*(n-1))
    r['MW_TOT'] = round(mw, 2)
    return r

