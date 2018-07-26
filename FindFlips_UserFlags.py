
# coding: utf-8

# In[21]:


__author__ = 'egatkinson'
#takes in 2 collapsed haploid bed files per individual, suffixed .A.bed and .B.bed, returns instances of double recombination events between them, which are assumed to be ancestry strand flip errors


# In[22]:


import argparse
import re
import numpy as np
import pprint


# In[27]:


def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--bed', help='path stem to collapsed bed files, usually ind name', required=True) #bedfile should be a list of bed files for indidivuals that have strands 1 and end ending with .A.bed and .B.bed
  parser.add_argument('--cM', help='cM window size for check', default=1) #default window size is start and stop between tracts is within 1cM

  args = parser.parse_args()
  return(args)                       


# In[29]:


#bedfile should be a list of bed files for indidivuals that have strands 1 and end ending with .A.bed and .B.bed
#it expects columns to be organized as in the collapse-ancestry.py script and will calculate switch errors according
#to cM position of the start and stop. Hard coded in to consider those within a 1cM window
args = parse_args()
cM = float(args.cM)
bedfile = open(args.bed)
beds = np.genfromtxt(bedfile, dtype=None)
## tallying up double recomb events where the start and stop is within 1cM of each other

for bed in beds:
    Atractsfile = open(bed + '.A.bed', 'r')
    Btractsfile = open(bed + '.B.bed', 'r')
    out = open(bed + '.DubRecomb.txt', 'w')
    Atracts = np.genfromtxt(Atractsfile, dtype=None)
    Btracts = np.genfromtxt(Btractsfile, dtype=None)

    for line in range(0, len(Atracts)):
        chromA = Atracts[line][0]
        startA = Atracts[line][5]
        endA = Atracts[line][7]
        ancA = Atracts[line][3]

        for item in range(0, len(Btracts)):
            chromB = Btracts[item][0]
            startB = Btracts[item][5]
            endB = Btracts[item][7]
            ancB = Btracts[item][3]

            count=0

            if chromA == chromB:
                if ancA != ancB:
                    #match to ensure that the ancestries are opposing. Right now, only modeling 2 way admixed scenario

                    #if the start and stop positions are within the window of each other consider this to be the same location
                    if (startA - cM) < startB < (startA + cM):
                        if (endA - cM) < endB < (endA + cM):
                            count = count +1
                        if count > 0:
                            output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chromA), str(startA), str(endA), str(ancA), str(startB), str(endB),str(ancB))
                            #output the chromosome number, start and end pos of A and B and their respective ancs for sanity check
                            #these are the double recombination events - the ones where there are matching tracts in opposing ancestries
                            out.write(output)
out.close()

