#!/usr/bin/env python
# Li Xue
# 12-Jul-2022 09:40

"""
Renumber model.pdb files based on ref.pdb (only keep the common part among all files)
NOTE: The input files are assumed to have match chain IDs (by `pdb_match_chn_batch.py`)

Usage:
        pdb_renum_batch.py <list of pdb files> <string:chain IDs>
Example:
        pdb_renum_batch.py example/file.list2 'M,N,O,P'

 INPUT (file.list: a file that contains a list of pdb files):
     path/ref.pdb
     path/model1.pdb
     path/model2.pdb

 NOTE about file.list:
   1. ref.pdb has to be on the first line of file.list
   2. path in file.list should be relative to the CWD or use absolute path


 OUTPUT: renumbered model.pdb files

"""

import os
from os.path import dirname
import sys
import subprocess
import re
from myfun import *

def check_input(args):

    if len(args) != 3:
        sys.stderr.write(__doc__)
        sys.exit(1)


def create_ini_commonPDB(refFL, chnIDs_ref, outDIR):
    # use ref.pdb as the initial common_A.pdb, common_B.pdb and so on
    # chnIDs_ref = ['A', 'B', 'C']
    for chnID in chnIDs_ref:
        outFL = f"{outDIR}/common_{chnID}.pdb"
        if os.path.isfile(outFL):
            os.remove(outFL)

        subprocess.check_call(f"cp {refFL} {outFL}", shell = True)

def rm_Xresidue(pdbFL):
    # if model.pdb has residues that are not in ref.pdb, pdb-pdbalign will use X in new_model.pdb to denote them
    # Delete residues that are not in ref.pdb, i.e., delete the ATOM lines with residue name of 'X'
    tmpFL = f"{pdbFL}.tmp"
    o = open(tmpFL,'w')
    f = open (pdbFL, 'r')
    for line in f:
        if (re.match('^ATOM', line)):
            if line[21] != "X":
                o.write(line)
    f.close()
    o.close()
    os.rename(tmpFL, pdbFL)
    print(f"X-residues are removed from {pdbFL}")

def create_commonPDB(refFL, pdbFLs, chnIDs_ref, outDIR):
    # compare pdbFLs with refFL and create common PDB files for each chain
    # the common PDB files are renumbered according to refFL

    create_ini_commonPDB(refFL, chnIDs_ref, outDIR)
    for chnID in chnIDs_ref:
        for pdbFL in pdbFLs:
            commonFL=f"{outDIR}/common_{chnID}.pdb"
            subprocess.check_call(f"pdb-pdbalign  {commonFL} {chnID} {pdbFL} {chnID} > {commonFL}_tmp", shell = True )
            # keep the new_model.pdb as the common_chnID.pdb
            subprocess.check_call(f"mv {commonFL}_tmp {commonFL}", shell = True)

def match_two_pdb(refPDB, chnID_ref, modelPDB, chnID_model, new_modelPDB):
    # match and renumber modelPDB according to refPDB
    print(f"match {modelPDB}_{chnID_model} to {refPDB}_{chnID_ref}...")
    subprocess.check_call(f"pdb-pdbalign {refPDB} {chnID_ref} {modelPDB} {chnID_model} > {new_modelPDB} ", shell = True )
    # clean up
    #subprocess.check_call(f"sed -i '/Warning/d' {new_pdbFL}", shell = True)
    rm_Xresidue(new_modelPDB)


check_input(sys.argv)

# read the list of pdb files
listFL = sys.argv[1]
chnIDs_ref = sys.argv[2].split(',')
pdbFLs = read_listFL(listFL)
refFL = pdbFLs.pop(0)
dir = dirname(pdbFLs[0])

# create common_A.pdb, common_B.pdb , common_chnID.pdb and so on
create_commonPDB(refFL, pdbFLs, chnIDs_ref, outDIR = dir)

# align and renumber each model.pdb to common_chnID.pdb
for pdbFL in pdbFLs:
    flname = basename(pdbFL, '.pdb')
    to_combine = []
    for chnID in chnIDs_ref:
        commonFL=f"{dir}/common_{chnID}.pdb"
        new_pdbFL = f"{dir}/{flname}_{chnID}_renum.pdb"
        match_two_pdb(commonFL, chnID, pdbFL, chnID, new_pdbFL)
        to_combine.append(new_pdbFL)

    # concatenate the chain files into one pdb file

    outFL = f"{dir}/{flname}_renum.pdb"
    #files = [ f"{dir}/{flname}_{chnID}_renum.pdb"  for chnID in chnIDs_ref]
    files = ' '.join(to_combine) # "file1 file2 file3"
    cat(files, outFL)
    print(f"{outFL} generated.\n")

    # clean up
    for i in to_combine: os.remove(i)
