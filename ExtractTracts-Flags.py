
# coding: utf-8

# In[100]:


#script to extract out ancestry segments from each reference population from admixed samples.
__author__ = 'egatkinson'


# In[101]:


import argparse
import re
import numpy as np


# In[ ]:


def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--msp', help='path stem to RFmix msp file, not including .msp.tsv', required=True)
  parser.add_argument('--vcf', help='path stem to RFmix input VCF with phased genotypes, not including .vcf suffix', required=True)
  args = parser.parse_args()
  return(args)    


# In[ ]:


args = parse_args()
mspfile = open(args.msp, + '.msp.tsv', 'r')
genofile = open(args.vcf + '.vcf', 'r')
out0 = open(args.genofile + '.anc0.vcf', 'w') #output for the extracted VCF anc 0
out1 = open(args.genofile + '.anc1.vcf', 'w') #output for the extracted VCF anc 1


# In[108]:


#save the first 6 lines, the VCF header
head = genofile.readlines()[:6]


# In[109]:


out0.write('\n'.strip().join(head))
out1.write('\n'.strip().join(head))
genofile = open(args.vcf + '.vcf', 'r') #have to re-open the VCF since this header logging closed it


# In[110]:


#initialize documenting the current window to check
chromosome = ("", 0, 0)


# In[111]:


for line in genofile:
    if line.startswith("#"):
        continue
    if not line:
        break #stop when get to the end of the file
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, genos = line.strip().split('\t', 9)
    genos = genos.replace('|', '\t').split('\t') #split each strand geno call apart from each other
    #genos = genos.split('\t').split('|')
    output0 = '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
    output1 = '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
    POS = int(POS)
    
    #optimized for quicker runtime - only move to next line when out of the current msp window
    #saves the current line until out of the window, then checks next line. VCF and window switches file should be in incremental order.
    while not (CHROM == chromosome[0] and ( chromosome[1] <= POS < chromosome[2])):
        ancs = mspfile.readline()
        if ancs.startswith("#"): #skip the header lines
            continue
        if not ancs:
            break #when get to the end of the msp file, stop
        chm, spos, epos, sgpos, egpos, nsnps, calls = ancs.strip().split('\t', 6)
        calls = calls.split('\t')
        chromosome = (chm, int(spos), int(epos))
        
    for i in range(len(genos)/2): #index by the number of individuals in the VCF file, should be the same number in the calls file
        genoA = str(genos[2*i])
        genoB = str(genos[2*i + 1])
        callA = str(calls[2*i])
        callB = str(calls[2*i + 1])

        #if the anc call is 0 (European CEU here), keep, replace 1 or other calls with missing data    
        if callA == '0':
            genoA0 = genoA
            genoA1 = "."
        elif callA == '1': 
            genoA0 = "."
            genoA1 = genoA   #if the anc call is 1 (African YRI here), keep, otherwise make into missing data

        if callB == '0':
            genoB0 = genoB
            genoB1 = "."
        elif callB == '1': 
            genoB0 = "."
            genoB1 = genoB   #if the anc call is 1 (African YRI here), keep, otherwise make into missing data

        output0 += '\t' + genoA0 + "|" + genoB0
        output1 += '\t' + genoA1 + "|" + genoB1
    output0 += '\n'
    output1 += '\n'
    out0.write(output0)
    out1.write(output1)
out0.close()
out1.close()

