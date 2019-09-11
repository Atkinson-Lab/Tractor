# coding: utf-8
#new version of flip catching/correcting script that detects strand flips and corrects them back in RFmix output
#this version fixes the RFmix output, .msp file with LAI calls for 2-way admixed individuals
__author__ = 'egatkinson'

import argparse
import re
import numpy as np

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--switches', help='file with LAI window calls regarding if genotypes need to be flipped there', required=True)
  parser.add_argument('--genofile', help='path stem to RFmix input VCF with phased genotypes, not including .vcf suffix', required=True)
  args = parser.parse_args()
  return(args)    

args = parse_args()
switches = open(args.switches, 'r')
genofile = open(args.genofile + '.vcf','r')
outVCF = open(args.genofile + '.unkinked.vcf', 'w')

#save the first 6 lines, the VCF header
#head = genofile.readlines()[:6]

#optimized to run much quicker
#outVCF.write('\n'.strip().join(head))
genofile = open(args.genofile + '.vcf')

#initialize documenting the current window to check
chromosome = ("", 0, 0)

#saves the current line in each file until out of the window, then checks next line. VCF and window switches file should be in incremental order.
for line in genofile:
    if line.startswith("#"):
        outVCF.write(line) 
        continue
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, genos = line.strip().split('\t', 9)
    genos = genos.split('\t')
    output = '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
    POS = int(POS)
    
    while not (CHROM == chromosome[0] and ( chromosome[1] <= POS < chromosome[2])):
        window = switches.readline()
        if not window:
            break #when get to the end of the file, stop
        chrom, startbp, stopbp, startcm, stopcm, snps, calls = window.strip().split('\t', 6)
        calls = calls.split('\t')
        chromosome = (chrom, int(startbp), int(stopbp))
        
    for i in range(len(genos)): #index by the number of individuals in the VCF file, should be the same number in the calls file
        hap = str(genos[i])
        newhap = hap

        if str(calls[i]) == "True":  #if the SNPs fall in a window which for that indiv is labeled switched, reverse the polymorphic sites
            if hap == "1|0":
                newhap = "0|1"
            elif hap == "0|1":
                newhap = "1|0"
        output += '\t' + newhap
    output += '\n'
    outVCF.write(output)
outVCF.close()          
