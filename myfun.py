import subprocess
import os
import re

def get_chain_list(pdbFL):
    # get a list of chain IDs from the input PDB pdbFL
    f = open(pdbFL, "r")
    chains0 = {}
    for line in f.readlines():
        if line[:4] != "ATOM": continue
        chains0[line[21]] = 1
        chains = []
    f.close()

    keys = chains0.keys()
    for key in keys:
        if key == " ": continue
        chains.append(key)
    return chains
def cat(files, outFL):
    # string: files = "file1 file2 file3"
    subprocess.check_call(f"cat {files} > {outFL}" , shell = True)

def read_listFL(listFL):
    f = open(listFL, 'r')
    pdbFLs = [line.rstrip() for line in f.readlines()]
    f.close()
    return pdbFLs

def basename(file, ext):
    return re.sub(ext, '', os.path.basename(file))
