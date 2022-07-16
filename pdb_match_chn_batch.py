#!/usr/bin/env python
"""
Map each chain ID in ref.pdb with a list of pdb files.
Minimum sequence identity: 70 %

Usage:
    pdb_match_chn_batch.py <a file of pdb file names>

Example:
    pdb_match_chn_batch.py example/file.list

INPUT (file.list: a file that contains a list of pdb files):
    path/ref.pdb
    path/model1.pdb
    path/model2.pdb
NOTE of file.list:
  1. ref.pdb has to be on the first line of file.list
  2. path in file.list should be relative to the CWD or use absolute path

Actions:
    1. obtain chain ID mappings between ref.pdb and to-be-mapped pdb files
    2. change the chain IDs of each pdb file to be mapped


OUTPUT:
    pdb files with matched chain IDs
"""

import sys, os
import pdb
import subprocess
from myfun import *

def check_input(args):

    if len(args) <1:
        sys.stderr.write(__doc__)
        sys.exit(1)

check_input(sys.argv[1:])
listFL = sys.argv[1]
pdbFLs = read_listFL(listFL)
refFL = pdbFLs.pop(0)

for pdbFL in pdbFLs:
    print (f"\n--- Mapping chain IDs of {pdbFL} to ref.pdb ---")

    # 1. generate chain ID map file
    dir = os.path.dirname(pdbFL)
    flname = basename(pdbFL,'.pdb')
    mapFL = f"{dir}/{flname}.chnMap"
    subprocess.check_call(f"pdb_match_chn.py {refFL} {pdbFL} > {mapFL}", shell=True)
    print (f"{mapFL} generated. Please double check!!!")

    #2. rename chain IDs of pdbFL
    outFL = f"{dir}/{flname}_newChnID.pdb"
    subprocess.check_call(f"pdb_rename_chain.py {pdbFL} {mapFL} {outFL}", shell=True)
    print (f"{outFL} with matched chain IDs generated.")



