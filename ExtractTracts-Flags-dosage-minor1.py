
# coding: utf-8

#script to extract out ancestry segments from each reference population from admixed samples.
__author__ = 'egatkinson'

import argparse
import re
import numpy as np


def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--msp', help='path stem to RFmix msp file, not including .msp.tsv', required=True)
  parser.add_argument('--vcf', help='path stem to RFmix input VCF with phased genotypes, not including .vcf suffix', required=True)
  args = parser.parse_args()
  return(args)    

args = parse_args()
mspfile = open(args.msp + '.msp.tsv', 'r')
genofile = open(args.vcf + '.vcf', 'r')
out0 = open(args.vcf + '.anc0.vcf', 'w') #output for the extracted VCF anc 0
out1 = open(args.vcf + '.anc1.vcf', 'w') #output for the extracted VCF anc 1
outdos0 = open(args.vcf + '.anc0.dosage.txt', 'w') #output dosages for each ancestry into separate files
outdos1 = open(args.vcf + '.anc1.dosage.txt', 'w') #output dosages for each ancestry into separate files

#save the first 6 lines, the VCF header
#head = genofile.readlines()[:6]
head = ""
for line in genofile:
	if line.startswith("#"):
		head = head + line
	else:
		break

out0.write('\n'.strip().join(head))
out1.write('\n'.strip().join(head))
genofile = open(args.vcf + '.vcf', 'r') #have to re-open the VCF since this header logging closed it


#initialize documenting the current window to check
chromosome = ("", 0, 0)


for line in genofile:
    if line.startswith("#"):
        continue
    if not line:
        break #stop when get to the end of the file
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, genos = line.strip().split('\t', 9)
    genos = genos.replace('|', '\t').split('\t') #split each strand geno call apart from each other
    #genos = genos.replace('/', '\t').split('\t') #split each strand geno call apart from each other, if didn't have phased pipes
    #genos = genos.split('\t').split('|')
    output0 = '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
    output1 = '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
    outputdos0 = '\t'.join([CHROM, POS, ID])
    outputdos1 = '\t'.join([CHROM, POS, ID])
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
        
    for i in range(len(genos)//2): #index by the number of individuals in the VCF file, should be the same number in the calls file
        genoA = str(genos[2*i])
        genoB = str(genos[2*i + 1])
        callA = str(calls[2*i])
        callB = str(calls[2*i + 1])
        count0 = 0
        count1 = 0
        count2 = 0

        #if the anc call is 0, keep, replace 1 or other calls with missing data    
        if callA == '0':
            genoA0 = genoA
            genoA1 = "."
            if genoA == '1':  #only tally up counts of the minor/risk allele for each ancestry; technically the alternate allele
            	count0 = count0 + 1
        elif callA == '1': 
            genoA0 = "."
            genoA1 = genoA   #if the anc call is 1, keep, otherwise make into missing data
            if genoA1 == '1':
            	count1 = count1 + 1

        if callB == '0':
            genoB0 = genoB
            genoB1 = "."
            if genoB0 == '1':
            	count0 = count0 + 1
        elif callB == '1': 
            genoB0 = "."
            genoB1 = genoB
            if genoB1 == '1':
            	count1 = count1 + 1

        output0 += '\t' + genoA0 + "|" + genoB0
        output1 += '\t' + genoA1 + "|" + genoB1
        outputdos0 += '\t' + str(count0)  #output the dosage for alt allele for each ancestry at each position for each indiv in the VCF file.
        outputdos1 += '\t' + str(count1)
    output0 += '\n'
    output1 += '\n'
    outputdos0 += '\n'
    outputdos1 += '\n'
    out0.write(output0)
    out1.write(output1)
    outdos0.write(outputdos0)
    outdos1.write(outputdos1)
out0.close()
out1.close()
outdos0.close()
outdos1.close()

