#!/usr/bin/env python
"""
Map each chain ID in ref.pdb with another pdb file.
The best match is kept. Therefore, one chain can be matched to multiple chains.
Minimum sequence identity: 70 %

Usage:
    pdb_match_chn.py <refe PDB file> <match PDB file>

Example:
    pdb_match_chn.py example/ref.pdb example/model1.pdb

OUTPUT:
    ref, modelpdb, score
    A, A, 98.56
    C, C, 100.0
    D, D, 90.86
    E, E, 84.01
"""

import sys, os
import pdb
import subprocess
import pdb
from myfun import *

def check_input(args):

    if not len(args):
        sys.stderr.write(__doc__)
        sys.exit(1)


check_input(sys.argv[1:])
min_identity = 70
pdbFL1 = sys.argv[1]
pdbFL2 = sys.argv[2]
chainsA = get_chain_list(pdbFL1) #chain IDs in pdbFL1
chainsB = get_chain_list(pdbFL2) #chain IDs in pdbFL2

for chainA in chainsA:
    match_nr = -1
    match_chain = "None"
    minscore = min_identity

    for chainBnr in range(0, len(chainsB)):
        chainB = chainsB[chainBnr]
        result = subprocess.check_output(f"pdb-pdbalignscore {pdbFL1} {chainA} {pdbFL2} {chainB}", shell=True)
        result = float(result)
        if result > minscore:
            minscore = result
            match_nr = chainBnr
            match_chain = chainsB[match_nr]
            matchScore = result
    print (f"{chainA}, {match_chain}, {matchScore}")
