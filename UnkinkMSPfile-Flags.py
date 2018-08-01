
# coding: utf-8

# In[ ]:


#new version of flip catching/correcting script that detects strand flips and corrects them back in RFmix output
#this version fixes the RFmix output, .msp file with LAI calls for 2-way admixed individuals
__author__ = 'egatkinson'


# In[67]:


import argparse
import re
import numpy as np


# In[ ]:


def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--msp', help='path stem to RFmix msp file, not including .msp.tsv', required=True)
  args = parser.parse_args()
  return(args)    


# In[ ]:


args = parse_args()
mspfile = open(args.msp + '.msp.tsv')
#beds = np.genfromtxt(bedfile, dtype=None)


# In[89]:


mspfile.readline()
header = mspfile.readline()
chrom, startbp, stopbp, startcm, stopcm, snps, samples = header.split('\t', 6)
state = {}
switched = {}

for i in range(len(samples)/2):
    state[i] = (float(0), [-1,-1])
    switched[i] = False


# In[90]:


out = open(args.msp + '.Unkinked.msp.tsv', 'w') #output will be input filename suffixed with "Unkinked"
out.write(header) #currently only 

for line in mspfile:
    chrom, startbp, stopbp, startcm, stopcm, snps, calls = line.strip().split('\t', 6)
    calls = calls.split('\t')
    output = '\t'.join([chrom, startbp, stopbp, startcm, stopcm, snps])
    for i in range(len(calls)/2):
        hapA = calls[2*i]
        hapB = calls[2*i + 1] 
        newhapA = hapA
        newhapB = hapB
        
        if hapA != hapB:            
            if [hapA,hapB] != state[i][1] and (float(startcm) - 1) < state[i][0]:
                switched[i] = not switched[i]
                
            if switched[i]:
                newhapA = hapB
                newhapB = hapA
           
            state[i] = (float(stopcm), [hapA,hapB])
        output += '\t' + newhapA + '\t' + newhapB
    output += '\n'
    out.write(output)
out.close()          

