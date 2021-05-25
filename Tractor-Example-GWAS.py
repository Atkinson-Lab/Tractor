#!/usr/bin/env python
# coding: utf-8

# # Initialize Hail
import argparse
import hail as hl
import numpy as np
hl.init()


#load in plotting features
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()


# # Load in the dosage files from Tractor
# ### Note: this will be the most time intensive step. Hail team is actively optimizing pieces of this infrastructure.
## user should modify the paths in the import steps to match the location (here shown for files on google cloud) of their datasets

#start loading in the ancestry 0 minor allele dosages
row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr, 'REF': hl.tstr, 'ALT': hl.tstr} 
anc0dos = hl.import_matrix_table('gs://.../Dataset.anc0.dosage.txt.gz', force_bgz=True, row_fields=row_fields, row_key=[], min_partitions=32) 
anc0dos = anc0dos.key_rows_by().drop('row_id')
anc0dos = anc0dos.key_rows_by(locus=hl.locus(anc0dos.CHROM, anc0dos.POS))     


#also load ancestry 1 allele dosages
row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr, 'REF': hl.tstr, 'ALT': hl.tstr} 
anc1dos = hl.import_matrix_table('gs://.../Dataset.anc1.dosage.txt.gz', force_bgz=True, row_fields=row_fields, row_key=[], min_partitions=32) 
anc1dos = anc1dos.key_rows_by().drop('row_id')
anc1dos = anc1dos.key_rows_by(locus=hl.locus(anc1dos.CHROM, anc1dos.POS))     

#Optional - save these temporary files to relieve  memory burden
anc0dos = anc0dos.checkpoint('gs://.../Dataset.anc0.dosage.mt')
anc1dos = anc1dos.checkpoint('gs://.../Dataset.anc1.dosage.mt')


#read in haplotype counts for the index ancestry, here anc0
row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr} 
hapcounts0 = hl.import_matrix_table('gs://.../Dataset.anc0.hapcount.txt.gz', force_bgz=True, row_fields=row_fields, row_key=[], min_partitions=32) 
hapcounts0 = hapcounts0.key_rows_by().drop('row_id')
hapcounts0 = hapcounts0.key_rows_by(locus=hl.locus(hapcounts0.CHROM, hapcounts0.POS))   

##if more than 2 way admixed, load in additional ancestry files with same sytax

#join the ancestry and haplotype dosage files together into a Hail matrix table
mt = anc0dos.annotate_entries(anc1dos = anc1dos[anc0dos.locus, anc0dos.col_id], hapcounts0 = hapcounts0[anc0dos.locus, anc0dos.col_id])

#Here's an example showing the implementation if cohort data is 3 way admixed
#mt = anc0dos.annotate_entries(anc1dos = anc1dos[anc0dos.locus, anc0dos.col_id], anc2dos = anc2dos[anc0dos.locus, anc0dos.col_id], hapcounts0 = hapcounts0[mt.locus, mt.col_id], hapcounts1 = hapcounts1[mt.locus, mt.col_id]))

## load in the phenotype phenotype and covariate file. The IDs should be consistent with those in the header of the Tractor dosage file
phenos = hl.import_table('gs://.../Dataset_pheno_and_covs.txt').key_by('IID')

#annotate the mt with the phenotypes
mt = mt.annotate_cols(phenos = phenos[mt.col_id])

#write out the newly annotated matrix table. Will run a lot faster if we load this in again after annotating things in due to Hail's lazy processing style.
mt.write('gs://.../Dataset.phenos_dosages.mt')


# # Run Tractor ancestry-aware linear regression

#read in the filtered and annotated mt file again to go faster
mt = hl.read_matrix_table('gs://.../Dataset.phenos_dosages.mt')

#Run GWAS on a single phenotype, here running on TC, Total Cholesterol, using the Tractor haplotype counts and dosages as described in the manuscript
#along with relevant covariates, here we show the covariates being sex, age, blood dilution level, and global ancestry fractions obtained from previous ADMIXTURE runs
#We are annotating the matrix table with the results of the GWAS. The results are saved as a struct here called TC
mt = mt.annotate_rows(TC = hl.agg.linreg(mt.TC,
                                                 [1.0, mt.hapcounts0.x, mt.anc0dos.x, mt.anc1dos.x, mt.isMale, mt.age, mt.dilutionFactor, hl.float(mt.admixFrac.AFR)]))


#to run phenotyps in batch, first make a list of the phenotype names. Here showing 4 blood lipid phenotypes
phenonames = ['TG', 'TC', 'HDLC', 'LDLC']

#run GWAS in batch for all listed phenotypes, annotating results in mt as separate structs named after the pheno
mt = mt.annotate_rows(results=hl.struct(**{pheno: hl.agg.linreg(mt[pheno], 
                                                                [1.0, mt.hapcounts0.x, mt.anc0dos.x, mt.anc1dos.x, mt.isMale, mt.age, mt.dilutionFactor, hl.float(mt.admixFrac.AFR)]) 
                                           for pheno in phenonames}))


#generate manhattan plots with bokeh, anc 0 and anc 1 are plotted seperately
#in the implementation above, the first (0 indexed) p value in the array is for the intercept, the second is ancestry haplotype counts at the site, the third is anc0, and the fourth is anc1, followed by the covariates
p = hl.plot.manhattan(mt.results.TC.p_value[2], title='Admixed Afr-Eur UKBB, TC, anc0', collect_all=False, significance_line=5e-08) 
show(p)

#plot up as a manhattan plot, anc 1
p = hl.plot.manhattan(mt.results.TC.p_value[3], title='Admixed Afr-Eur UKBB, TC, anc1', collect_all=False, significance_line=5e-08) 
show(p)


#make a QQ plot for TC anc0
p = hl.plot.qq(mt.results.TC.p_value[2], title="QQ plot, TC, anc0")
show(p)

#make a QQ plot for TC anc1
p = hl.plot.qq(mt.results.TC.p_value[3], title="QQ plot, TC, anc1")
show(p)

