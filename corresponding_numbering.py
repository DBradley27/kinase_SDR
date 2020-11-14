#!/usr/bin/python

import os
import sys

position_list = []
file_in = open(sys.argv[1],"r")
file_out = open(sys.argv[2],"w")

while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    position_list = line.split(", ")
    break
file_in.close()


#for i in range(len(position_list)):
    #print i+1, int(position_list[i])+1

for i in range(len(position_list)):
    a = i+1
    b = int(position_list[i])+1
    c = str(a)+" "+str(b)
    file_out.write(c+"\n")

file_out.close()
