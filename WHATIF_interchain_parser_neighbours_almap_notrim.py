#!/usr/bin/python

# This is a script that maps the position a site in the trimmed alignment to the position in the original, untrimmed alignment. 

import os
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

# Input files.

Neighbours_in = open(sys.argv[1],"r")
alignment_mapping_in = open(sys.argv[2],"r")
file_out = open(sys.argv[5],"w")

CDK2_neighbours = pd.read_csv(Neighbours_in,sep=" ",header=None)
mapping = pd.read_csv(alignment_mapping_in,sep=" ",header=None)

CDK2_neighbours.columns = ['column_one','column_two']
CDK2_neighbours = CDK2_neighbours.dropna(subset = ['column_one', 'column_two'])

# Load sequence alignment of protein domains; create dictionary (sequence_dict) mapping protein name to protein sequence

sequence_dict = {}
for record in SeqIO.parse(open(sys.argv[3]), "fasta"):
    gene_name = record.id.split("-")[0]
    sequence = record.seq.tostring().replace(" ","")
    sequence_dict[gene_name] = sequence
    if "U" in sequence:
        print tab[0]
        sys.exit()

# Generate dictionary (corresponding table) for each protein that maps position in alignment to position in the sequence
#  (N.B. position relative to start of domain, not raw sequence).

# Assign each dictionary to the protein name in a higher-level dictionary (corresponding_table_dict).

corresponding_table_dict = {}

for gene in sequence_dict:
    i, j = -1, -1
    corresponding_table = {}
    for aa in sequence_dict[gene]:
        i=i+1
        if aa != "-":
            j = j+1
            corresponding_table[i] = j
    corresponding_table_dict[gene] = corresponding_table

# Adjust the alignment mapping to account for the starting position of the domain in the raw sequence.

for i in corresponding_table_dict:
    start = 1    
    for j in corresponding_table_dict[i]:         
        corresponding_table_dict[i][j] = corresponding_table_dict[i][j] + start 

# Output mapping

text = str(corresponding_table_dict)
gene_id = text[text.find(sys.argv[4]):text.find(sys.argv[4])+21]
gene_id = gene_id[:gene_id.find("'")]

print gene_id


aln_mapping = corresponding_table_dict[gene_id]
aln_mapping = pd.Series(aln_mapping)


aln_mapping_fin = {}

for i in mapping.iloc[:,1]:
      if math.isnan(float(i)) == False:
          i = int(i)
      s = aln_mapping.get(i-1,'Gap')
      aln_mapping_fin[i] = s  
    
aln_mapping_fin = pd.DataFrame(aln_mapping_fin.items()) 
aln_mapping_fin_sort = aln_mapping_fin.sort([0],ascending=True)
aln_mapping_fin_sort.insert(0,'aln',list(mapping.iloc[:,0]))

aln_mapping_fin_sort.columns = ['aln','trim_map','seq_map']



row = aln_mapping_fin_sort[aln_mapping_fin_sort['seq_map'].isin([1000])]


for i in range(len(CDK2_neighbours)):
    row = aln_mapping_fin_sort[aln_mapping_fin_sort['seq_map'].isin([int(CDK2_neighbours.iloc[i,0])])]
    if len(row) == 1:
        number = int(row.iloc[:,0])
        CDK2_neighbours.iloc[i,0] = number
    else:
        CDK2_neighbours.iloc[i,0] = 9999 

for i in range(len(CDK2_neighbours)):
    row = aln_mapping_fin_sort[aln_mapping_fin_sort['seq_map'].isin([int(CDK2_neighbours.iloc[i,1])])]
    if len(row) == 1:
        number = int(row.iloc[:,0])
        CDK2_neighbours.iloc[i,1] = number
    else:
        CDK2_neighbours.iloc[i,1] = 9999


CDK2_neighbours = CDK2_neighbours.astype(int)
CDK2_neighbours.to_csv(sys.argv[5],sep=" ",header=None, index=False)





Neighbours_in.close()
alignment_mapping_in.close()
file_out.close()



