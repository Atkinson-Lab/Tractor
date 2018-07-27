
# coding: utf-8

# # A script for flipping strand flips induced by phasing, in (mostly) collapsed bed files

# In[1]:


__author__ = 'egatkinson'


# In[2]:


import argparse
import numpy as np
import re
import pprint


# In[6]:


def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--bed', help='path stem to collapsed bed files, usually ind name', required=True) #bedfile should be a list of bed files for indidivuals that have strands 1 and end ending with .A.bed and .B.bed
  args = parser.parse_args()
  return(args)    


# In[7]:


#bedfile should be a list of bed files for indidivuals that have strands 1 and end ending with .A.bed and .B.bed
#it expects columns to be organized as in the collapse-ancestry.py script and will calculate switch errors according to cM position of the start and stop.
#output will be bed line from input suffixed with .flipped.*.bed
args = parse_args()
bedfile = open(args.bed)
beds = np.genfromtxt(bedfile, dtype=None)
# IF the entry matches one in the double recomb list, then print the opposite ancestry for that tract 
for bed in beds:
    Atractsfile = open(bed + '.A.bed', 'r')
    Btractsfile = open(bed + '.B.bed', 'r')
    dubrecomfile = open(bed + '.DubRecomb.txt', 'r') #assumes the double recomb file was produced from the previous script
    outA = open(bed + '.flipped.A.bed', 'w')
    outB = open(bed + '.flipped.B.bed', 'w')
    Atracts = np.genfromtxt(Atractsfile, dtype=None)
    Btracts = np.genfromtxt(Btractsfile, dtype=None)
    dubrecom = np.genfromtxt(dubrecomfile, dtype=None)
    
    for line in range(0, len(Atracts)):
        for entry in range(0, len(dubrecom)):
            if dubrecom[entry][0] == Atracts[line][0]: #if the chromosomes are the same
                if dubrecom[entry][1] == Atracts[line][5]:
                    Atracts[line][3] = dubrecom[entry][6]
        outputA = '\t'.join(map(str, Atracts[line])) + '\n'
                    #if the A file cM start position matches one of the logged double recomb even sites, then print the oppsite anc from the other strand
        outA.write(outputA)
    
    for line in range(0, len(Btracts)): #And do again for the B bed file
        for entry in range(0, len(dubrecom)):
            if dubrecom[entry][0] == Btracts[line][0]: #if the chromosomes are the same
                if dubrecom[entry][4] == Btracts[line][5]:
                    Btracts[line][3] = dubrecom[entry][3]
        outputB = '\t'.join(map(str, Btracts[line])) + '\n'
        outB.write(outputB)
    

outB.close()
outA.close()

