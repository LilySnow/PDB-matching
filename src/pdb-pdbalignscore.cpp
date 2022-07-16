//Aligns PDB 2 to PDB 1 and renumbers accordingly, provided that all residues of PDB 2 are present in PDB 1. If they are not present, a warning is generated, and the residue is named X 1, 2, ...

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <deque>
#include <map>
#include <algorithm>
using namespace std;

inline char getresidue(char code[3]) {
    if (!strcmp(code, "ALA")) return 'A';
    if (!strcmp(code, "CYS")) return 'C';
    if (!strcmp(code, "ASP")) return 'D';
    if (!strcmp(code, "GLU")) return 'E';
    if (!strcmp(code, "PHE")) return 'F';
    if (!strcmp(code, "GLY")) return 'G';
    if (!strcmp(code, "HIS")) return 'H';
    if (!strcmp(code, "ILE")) return 'I';
    if (!strcmp(code, "LYS")) return 'K';
    if (!strcmp(code, "LEU")) return 'L';
    if (!strcmp(code, "MET")) return 'M';
    if (!strcmp(code, "ASN")) return 'N';
    if (!strcmp(code, "PRO")) return 'P';
    if (!strcmp(code, "GLN")) return 'Q';
    if (!strcmp(code, "ARG")) return 'R';
    if (!strcmp(code, "SER")) return 'S';
    if (!strcmp(code, "THR")) return 'T';
    if (!strcmp(code, "VAL")) return 'V';
    if (!strcmp(code, "TRP")) return 'W';
    if (!strcmp(code, "TYR")) return 'Y';
    return 'X';
}

int extractpdb(const char *pdb, char chainid, char *&seq, int *&num) {
  bool has_chain = 0;
  FILE *fil = fopen(pdb, "r");
  if (fil == NULL) {fprintf(stderr, "PDB file %s does not exist\n", pdb); exit(1); }
  int num0[10000];
  char seq0[10000];
  int counter = 0;
  int curr_resnr = -1000;

  while (!feof(fil)) {
    char buf[100];
    char code[6];
    char name[3];
    fgets(buf, 100, fil);
    sscanf(buf, "%s %*s %*s %s", code, name);
    if (!strncmp(code,"ATOM", 4)) {
      char currchainid = buf[21];
      if (currchainid != chainid) continue;
      has_chain = 1;
      int resnr = atoi(buf+22);
      if (resnr != curr_resnr) {
        seq0[counter] = getresidue(name);
	num0[counter] = resnr;
	curr_resnr = resnr;
	counter++;
      }
    }
  }
  fclose(fil);
  if (!counter) {fprintf(stderr, "%s is not a PDB file\n", pdb); exit(1); }
  if (!has_chain) {fprintf(stderr, "%s has no chainid %c\n", chainid); exit(1);}
  seq = new char[counter+1];
  memcpy(seq, seq0, counter);
  seq[counter] = 0;
  num = new int[counter];
  memcpy(num, num0, counter *sizeof(int));
  return counter;
}

int match = 40;
int mismtch = -20;
int gap_open = -100;
int gap_cont = -2;


struct AlignmentPos {
  int score;
  int direc; //0 = left, 1 = up, 2 = left-up
};

void initialize(AlignmentPos **aln, int seqlen1, int seqlen2){
  int n;
  for (n = 0; n < seqlen1 + 1; n++) {
    for (int nn = 0; nn < seqlen2 + 1; nn++) aln[nn][n].score = 0;
  }
  aln[0][0].score = 0;
  aln[0][0].direc = 0;
}

bool feedforward(AlignmentPos **aln, int seqlen1, int seqlen2, int row) {
  int n;
  bool change = 0;
  for (n = 0; n < seqlen1; n++) {
    int difscore = gap_open;
    if (n > 0 && aln[row][n].direc == 0) difscore = gap_cont;
    int newscore = aln[row][n].score + difscore;
    if (newscore > aln[row][n+1].score) {
      change = 1;
      aln[row][n+1].score = newscore;
      aln[row][n+1].direc = 0;
    } 
  }
  return change;
}

bool feeddown(AlignmentPos **aln, int seqlen1, int seqlen2,int column) {
  int n;
  bool change = 0;
  for (n = 0; n < seqlen2; n++) {
    int difscore = gap_open;
    if (n > 0 && aln[n][column].direc == 1) difscore = gap_cont;
    int newscore = aln[n][column].score + difscore;
    if (newscore > aln[n+1][column].score) {
      change = 1;
      aln[n+1][column].score = newscore;
      aln[n+1][column].direc = 1;
    } 
  }
  return change;
}
bool feedforwarddown(AlignmentPos **aln, char *seq1, char *seq2, int seqlen1, int seqlen2, bool startleft, int diagonal) {
  bool change = 0;
  int start1 = 0, start2 = 0;
  if (startleft) start1 = diagonal; else start2 = diagonal;
  for (int n1 = start1, n2 = start2; n1 < seqlen2 && n2 < seqlen1; n1++, n2++) {
    int difscore = mismtch;
    if (seq2[n1] == seq1[n2]) difscore = match;
    int newscore = aln[n1][n2].score + difscore;
    if (newscore > aln[n1+1][n2+1].score) {
      change = 1;
      aln[n1+1][n2+1].score = newscore;
      aln[n1+1][n2+1].direc = 2;
    }
  }
  return change;
}

int main(int argc, char *argv[]) {
  int n;

  if (argc < 5) {
    cerr << "Usage: pdb-pdbalignscore <pdb file 1> <chainid 1> <pdb file 2> <chainid2> " << endl;
    return 1;
  }
  char *seq1; int *num1;
  char *seq2; int *num2;
  
  char chainid1 = argv[2][0];
  char chainid2 = argv[4][0];
  int seqlen1 = extractpdb(argv[1], chainid1, seq1, num1);
  int seqlen2 = extractpdb(argv[3], chainid2, seq2, num2);
  AlignmentPos **aln = new AlignmentPos*[seqlen2+1]; //seqlen2+1 rows, ranging from 0 to seqlen
  for (n = 0; n < seqlen2 + 1; n++) {
    aln[n] = new AlignmentPos[seqlen1+1];
  }
  
  initialize(aln, seqlen1, seqlen2);
  bool change = 1;
  int counter = 0;
  while (change) {
    change = 0;
    for (n = 0; n < seqlen2; n++) {
      if (feedforwarddown(aln, seq1, seq2, seqlen1, seqlen2, 1, n)) change = 1;
    }
    for (n = 1; n < seqlen1; n++) {
      if (feedforwarddown(aln, seq1, seq2, seqlen1, seqlen2, 0, n)) change = 1;
    }    
    for (n = 0; n < seqlen2+1; n++) {
      if (feedforward(aln,seqlen1, seqlen2, n)) change = 1;
    }    
    for (n = 0; n < seqlen1+1; n++) {
      if (feeddown(aln,seqlen1, seqlen2, n)) change = 1;
    }        
    counter++;
  }
  double score = double(100) * aln[seqlen2][seqlen1].score;
  int seqlenmin = seqlen1; if (seqlen2 < seqlen1) seqlenmin = seqlen2;
  score = score / match / seqlenmin;
  printf("%.2f\n", score);
  return 0;
}
