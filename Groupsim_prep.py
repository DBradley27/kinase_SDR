#!/usr/bin/python

# This is a script to make the full-length sequence identifiers from UniProt Fasta files more readable.
# It concatenates three elements from the original identifier: the accession, the entry name, and the protein name
# I  this instance I have removed reference to the species origin "e.g. HUMAN", to
# eliminate redundancy assuming that the user is working with sequences from a single species.

import os
import sys

# Two arguments: read file and write file. Assign command line arguments as variables in Python.

file_in = open(sys.argv[1],"r")
file_out = open(sys.argv[3],"w")

# While loop: runs consecutively through each line in read file, breaks when the end of the document is reached

while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    if line[0] == ">":

        line = line.replace("(","")
        line = line.replace(")","")
        tab = line.split("/")
        name = tab[0]
        
        record = name+"|"+sys.argv[2]

        line = record
#+"/"+tab[1]

    file_out.write(line+"\n")


file_in.close()
file_out.close()
