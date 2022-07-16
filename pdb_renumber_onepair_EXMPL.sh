# This script maps and renumbers one model.pdb according to ref.pdb
# OUTPUT:
#    a new model.pdb with renumbering
# Note:
# The input pdb files have to have seg id
#
# Here is an example:  Suppose we want to match ref.pdb and chain M with model.pdb with chain A

echo "1. Match and renumber model.pdb to ref.pdb"
chain_ID_ref='M'
chain_ID_mobile = 'A'
pdb-pdbalign ref.pdb $chain_ID_ref model.pdb $chain_ID_mobile > model_new.pdb #if model.pdb has residues that are not in ref.pdb, model_new.pdb will use X to denote them

# delete the warning line in common.pdb
sed -i '/Warning/d' model_new.pdb

# optional: Delete residues that are not in ref.pdb, i.e., delete the ATOM lines with residue name of 'X'
echo "2. (optional) Delete residues that are in model.pdb but not in ref.pdb"
awk '$1 == "ATOM" && substr($0,22,1) != "X"' model_new.pdb  >  common_${chain_ID_ref}.pdb

echo "model_new.pdb generated"
