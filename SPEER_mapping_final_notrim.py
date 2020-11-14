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

SPEER_in = open(sys.argv[1],"r")
alignment_mapping_in = open(sys.argv[2],"r")
file_out = open(sys.argv[5],"w")



# Mapping from trimmed alignment to full alignment; incorporate in SPEER data frame #22 for apirwise ## for 8

SPEER = pd.read_table(SPEER_in,delim_whitespace=True,skiprows=22,header=None)

#SPEER = SPEER.drop(SPEER.columns[0], axis=1)
print SPEER
#SPEER = SPEER.drop(SPEER.columns[5], axis=1)
#SPEER = SPEER.drop(SPEER.columns[5], axis=1)
#SPEER = SPEER.drop(SPEER.columns[5], axis=1)

mapping = pd.read_csv(alignment_mapping_in,sep=" ",header=None)

SPEER.columns = ['Column','%ID','SPEER-score','Z-score','P-value']

Columns = list(SPEER.iloc[:,0])

SPEER.iloc[:,0] = [int(x)+1 for x in SPEER.iloc[:,0]]

SPEER = SPEER.sort('Column',ascending=True)

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
gene_id = text[text.find(sys.argv[4]):text.find(sys.argv[4])+30]
print gene_id
print gene_id
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

#for i in range(len(column_three)):
    #if column_three[i] != "Gap":
        #column_three[i] = int(column_three[i]) - 1

#print SPEER.iloc[:,0]

SPEER_col = map(int, list(SPEER.iloc[:,0]))

crossover = set(list(aln_mapping_fin_sort.iloc[:,0])).intersection(SPEER_col)

crossover = list(crossover)

aln_mapping_truncated = pd.DataFrame()

for i in range(len(crossover)):
    for j in range(len(aln_mapping_fin_sort)):
    	if crossover[i] == aln_mapping_fin_sort.iloc[j,0]:
        	aln_mapping_truncated = aln_mapping_truncated.append(aln_mapping_fin_sort.iloc[j,:],ignore_index=True)    
	    	    



column_two = pd.DataFrame(aln_mapping_truncated.iloc[:,0])
column_three = pd.DataFrame(aln_mapping_truncated.iloc[:,1])

#column_two = column_two.sort([0],ascending=True)

#column_three = column_three.sort([1],ascending=True)
#print column_two.shape
#print column_three.shape
#print SPEER.shape

SPEER.insert(1, 'trim_map',list(column_two.iloc[:,0]))

SPEER.insert(2,'aln_map_'+sys.argv[4], list(column_three.iloc[:,0]))

SPEER.to_csv(sys.argv[5], sep="\t")


SPEER_in.close()
alignment_mapping_in.close()
file_out.close()



