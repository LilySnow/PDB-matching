# PDB-matching


Before we calculate RMSD, it is often a good pratice to make sure the ref.pdb
and model.pdb files have the same chain IDs and the same residue numbering.

## Prerequisite
1. install pdb-tools by HADDOCK: https://github.com/haddocking/pdb-tools
2. Add this line to `~/.bashrc`:
   ```
       export PATH="$PATH:YOUR_DIR/PDB-matching"
   ```
   Replace YOUR_DIR with your local folder for PDB-matching.
   Then in your terminal run:
   ```source ~/.bashrc```
   

## Two main scripts
- Step 1. `pdb_match_chn_batch.py`: match chain IDs of pdb files to ref.pdb.
  Output `_newChnID.pdb` files.

  *Note: This step can be skipped if model.pdb files have already matched chain IDs.*

- Step 2. `pdb_renum_batch.py`: align and renumber pdb files to ref.pdb. Output
  `_renum.pdb` files.

  *Note: Update `file.list` to `_newChnID.pdb` files.*

## Step 1: Match chain IDs to ref.pdb

### For two PDB files
- `pdb_match_chn.py`:

    -> match chain IDs for two PDB files (Does not always work well. Need to
    manually check the results!)
    
- `pdb_rename_chain.py`:

    -> rename chain IDs of a pdb file according to the ref.pdb (the chain ID map file generated by `pdb_match_chn.py`)

### For a list of pdb files
- `pdb_match_chn_batch.py`


## Step 2: Renumber pdb files to ref.pdb

### For two PDB files
- `pdb-pdbalign ref.pdb chainID model.pdb chainID`:  

    -> align and renumber model.pdb according to ref.pdb. Output: new model.pdb
    
    -> If model.pdb has a residue that is not in ref.pdb, new_model.pdb will use X to denote this residue: e.g., `ATOM   1447  N   GLN X   1 ...`.
- `pdb_renumber_onepair_EXMPL.sh`:  

    -> an example code

### For a list of pdb files
- `pdb_renum_batch.py file.list 'A,C,D,E'`, where 'A,C,D,E' are matched chain IDs

