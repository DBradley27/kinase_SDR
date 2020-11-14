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

GroupSim_in = open(sys.argv[1],"r")
alignment_mapping_in = open(sys.argv[2],"r")
file_out = open(sys.argv[5],"w")



# Mapping from trimmed alignment to full alignment; incorporate in GroupSim data frame

GroupSim = pd.read_csv(GroupSim_in,sep="\t",skiprows=4,header=None)

mapping = pd.read_csv(alignment_mapping_in,sep=" ",header=None)
trim_map = mapping.iloc[:,1]
trim_map = [x - 0 for x in trim_map]
trim_map = pd.DataFrame(trim_map)

GroupSim.insert(1,'trim_map', trim_map)

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

column_three = aln_mapping_fin_sort.iloc[:,1]
column_three = list(column_three)

#for i in range(len(column_three)):
    #if column_three[i] != "Gap":
        #column_three[i] = int(column_three[i]) - 1

column_three = pd.DataFrame(column_three)

print column_three

GroupSim.insert(2,'aln_map_'+sys.argv[4],column_three)

GroupSim = GroupSim.iloc[:,0:4]

GroupSim.to_csv(sys.argv[5],sep="\t")


GroupSim_in.close()
alignment_mapping_in.close()
file_out.close()



