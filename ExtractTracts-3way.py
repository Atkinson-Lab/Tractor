# coding: utf-8
#Tractor script to extract out ancestry segments from each reference population from admixed samples.
__author__ = 'egatkinson'


import argparse
import re
import numpy as np

USAGE = """
ExtractTracts-Flags-multi.py --msp  <an ancestral calls file produced by RFmix version 2, suffixed with .msp.tsv>
							 --vcf <VCF file suffixed with .vcf>
"""

#input is expected to be a VCF file suffixed with .vcf and an ancestral calls file produced by RFmix version 2, suffixed with .msp.tsv
#output will be returned in the same location as input VCF file with same naming convention

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
out2 = open(args.vcf + '.anc2.vcf', 'w') #output for the extracted VCF anc 2
outdos0 = open(args.vcf + '.anc0.dosage.txt', 'w') #output dosages for each ancestry into separate files
outdos1 = open(args.vcf + '.anc1.dosage.txt', 'w') #output dosages for each ancestry into separate files
outdos2 = open(args.vcf + '.anc2.dosage.txt', 'w') #output dosages for each ancestry into separate files
outancdos0 = open(args.vcf + '.anc0.hapcount.txt', 'w')  # output number of haplotype for each ancestry into separate files
outancdos1 = open(args.vcf + '.anc1.hapcount.txt', 'w')  # output number of haplotype for each ancestry into separate files
outancdos2 = open(args.vcf + '.anc2.hapcount.txt', 'w')  # output number of haplotype for each ancestry into separate files

#save the VCF header
head = ""
for line in genofile:
	if line.startswith("#"):
		head = head + line
	else:
		break

out0.write('\n'.strip().join(head))
out1.write('\n'.strip().join(head))
out2.write('\n'.strip().join(head))
genofile = open(args.vcf + '.vcf', 'r') #may have to re-open the VCF since this header logging closed it


#initialize documenting the current window to check
chromosome = ("", 0, 0)

for line in genofile:
	if line.startswith("#"): #save the header of the VCF file
		continue
	if not line:
		break #stop when get to the end of the file
	CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, genos = line.strip().split('\t', 9)
	genos = genos.replace('|', '\t').split('\t') #split each strand geno call apart from each other
	#genos = genos.split('\t').split('|')
	output0 = '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
	output1 = '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
	output2 = '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT])
	outputdos0 = '\t'.join([CHROM, POS, ID, REF, ALT])
	outputdos1 = '\t'.join([CHROM, POS, ID, REF, ALT])
	outputdos2 = '\t'.join([CHROM, POS, ID, REF, ALT])
	outputancdos0 = '\t'.join([CHROM, POS, ID, REF, ALT])
	outputancdos1 = '\t'.join([CHROM, POS, ID, REF, ALT])
	outputancdos2 = '\t'.join([CHROM, POS, ID, REF, ALT])
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
		count0 = 0
		count1 = 0
		count2 = 0
		count_anc0: int = 0
		count_anc1: int = 0

		#if the anc call is 0, keep, replace 1 or other calls with missing data
		if callA == '0':
			genoA0 = genoA
			genoA1 = "."
			genoA2 = "."
			count_anc0 = count_anc0 + 1  # tally up the ancestral haplotypes present at each site
			if genoA0 == '1':
				count0 = count0 + 1
		elif callA == '1':
			genoA0 = "."
			genoA1 = genoA   #if the anc call is 1, keep, otherwise make into missing data
			genoA2 = "."
			count_anc1 = count_anc1 + 1
			if genoA1 == '1':
				count1 = count1 + 1
		elif callA == '2': #and also get the dosage and pieces for anc call 2
			genoA0 = "."
			genoA1 = "."
			genoA2 = genoA
			count_anc2 = count_anc2 + 1
			if genoA2 == '1':
				count2 = count2 + 1

		if callB == '0':
			genoB0 = genoB
			genoB1 = "."
			genoB2 = "."
			count_anc0 = count_anc0 + 1
			if genoB0 == '1':
				count0 = count0 + 1
		elif callB == '1':
			genoB0 = "."
			genoB1 = genoB
			genoB2 = "."
			count_anc1 = count_anc1 + 1
			if genoB1 == '1':
				count1 = count1 + 1
		elif callB == '2':
			genoB0 = "."
			genoB1 = "."
			genoB2 = genoB
			count_anc2 = count_anc2 + 1
			if genoB2 == '1':
				count2 = count2 + 1

		output0 += '\t' + genoA0 + "|" + genoB0
		output1 += '\t' + genoA1 + "|" + genoB1
		output2 += '\t' + genoA2 + "|" + genoB2
		outputdos0 += '\t' + str(count0)
		outputdos1 += '\t' + str(count1)
		outputdos2 += '\t' + str(count2)
		outputancdos0 += '\t' + str(count_anc0)
		outputancdos1 += '\t' + str(count_anc1)
		outputancdos2 += '\t' + str(count_anc2)

	output0 += '\n'
	output1 += '\n'
	output2 += '\n'
	outputdos0 += '\n'
	outputdos1 += '\n'
	outputdos2 += '\n'
	outputancdos0 += '\n'
	outputancdos1 += '\n'
	outputancdos2 += '\n'
	out0.write(output0)
	out1.write(output1)
	out2.write(output2)
	outdos0.write(outputdos0)
	outdos1.write(outputdos1)
	outdos2.write(outputdos2)
	outancdos0.write(outputancdos0)
	outancdos1.write(outputancdos1)
	outancdos2.write(outputancdos2)
out0.close()
out1.close()
out2.close()
outdos0.close()
outdos1.close()
outdos2.close()
outancdos0.close()
outancdos1.close()
outancdos2.close()