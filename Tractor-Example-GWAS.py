#!/usr/bin/env python
# coding: utf-8

# # Initialize Hail

# In[1]:


import argparse
import hail as hl
import numpy as np
hl.init()


# In[2]:


#load in plotting features
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()


# # Load in the filtered and phenotype annotated genotype data for the individuals
# ### Note: the example file is in Matrix Table format, the native Hail format. VCF formats may also be imported with 'hl.import_vcf'. See documentation for more details: https://hail.is/

# In[ ]:


#key columns by sample and rows by rsID for easier merge with dosage data
mt = hl.read_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/Results/UKBB_AfEur_QCed_lipids2.mt').key_rows_by('locus')


# # Load in the dosage files from Tractor
# ### Note: this will be the most time intensive step. Hail team is actively optimizing pieces of this infrastructure.

# In[ ]:


row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr, 'REF': hl.tstr, 'ALT': hl.tstr} 
anc0dos = hl.import_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.autosomes.anc0.dosage_v1.txt.gz', 
                                 force_bgz=True, row_fields=row_fields, row_key=[], min_partitions=32) 
anc0dos = anc0dos.key_rows_by().drop('row_id')
anc0dos = anc0dos.key_rows_by(locus=hl.locus(anc0dos.CHROM, anc0dos.POS))     


# In[ ]:


row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr, 'REF': hl.tstr, 'ALT': hl.tstr} 
anc1dos = hl.import_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.autosomes.anc1.dosage_v1.txt.gz', 
                                 force_bgz=True, row_fields=row_fields, row_key=[], min_partitions=32) 
anc1dos = anc1dos.key_rows_by().drop('row_id')
anc1dos = anc1dos.key_rows_by(locus=hl.locus(anc1dos.CHROM, anc1dos.POS))     


# In[ ]:


#Optional - save these temporary files to relieve  memory burden
anc1dos = anc1dos.checkpoint('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.autosomes.anc1.dosage_v1.mt')
anc0dos = anc0dos.checkpoint('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.autosomes.anc0.dosage_v1.mt')


# In[ ]:


#read in haplotype counts for anc0, here African, the index ancestry
row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr} 
hapcounts0 = hl.import_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.auto.anc0.hapcount.txt.gz', 
                                 force_bgz=True, row_fields=row_fields, row_key=[], min_partitions=32) 
hapcounts0 = hapcounts0.key_rows_by().drop('row_id')
hapcounts0 = hapcounts0.key_rows_by(locus=hl.locus(hapcounts0.CHROM, hapcounts0.POS))   


# In[ ]:


#join the dosage files to the genotype data. 
#Specifically, annotating the samples with their info for how many copies of the minor allele per ancestry were seen and index haplotype counts
mt = mt.annotate_entries(anc0dos = anc0dos[mt.locus, mt.s], anc1dos = anc1dos[mt.locus, mt.s], hapcounts0 = hapcounts0[mt.locus, mt.s])


# In[ ]:


#write out the newly annotated matrix table. Will run a lot faster if we load this in again after annotating things in due to Hail processing style.
mt.write('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids_dosages.mt')


# # Run linear regression across all phenotypes.

# In[3]:


#read in the filtered and annotated mt file again to go faster
mt = hl.read_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids_dosages.mt')


# In[4]:


#Run GWAS on a single phenotype, here running on TC, Total Cholesterol, using the Tractor haplotype counts and dosages
#along with the covariates of ages, sex, blood dilution level, and global ancestry fractions obtained from previous ADMIXTURE runs
mt = mt.annotate_rows(TC = hl.agg.linreg(mt.TC,
                                                 [1.0, mt.hapcounts0.x, mt.anc0dos.x, mt.anc1dos.x, mt.isMale, mt.age, mt.dilutionFactor, hl.float(mt.admixFrac.AFR)]))


# In[45]:


#to run phenotyps in batch, make a list of the phenotype names
phenonames = ['TG', 'TC', 'HDLC', 'LDLC']


# In[5]:


#run GWAS in batch for all listed phenotypes, annotating with the pheno name
mt = mt.annotate_rows(results=hl.struct(**{pheno: hl.agg.linreg(mt[pheno], 
                                                                [1.0, mt.hapcounts0.x, mt.anc0dos.x, mt.anc1dos.x, mt.isMale, mt.age, mt.dilutionFactor, hl.float(mt.admixFrac.AFR)]) 
                                           for pheno in phenonames}))


# In[ ]:


#plot up as a manhattan plot, anc 0 and anc 1 seperately
#in the implementation above, the first p value in the array is for the intercept, the second is ancestry counts at the site, the third is anc0, and the fourth is anc1, followed by the covariates
p = hl.plot.manhattan(mt.results.TC.p_value[2], title='Admixed Afr-Eur UKBB, TC, anc0', collect_all=False, significance_line=5e-08) 
show(p)


# In[ ]:


#plot up as a manhattan plot, anc 1
p = hl.plot.manhattan(mt.results.TC.p_value[3], title='Admixed Afr-Eur UKBB, TC, anc1', collect_all=False, significance_line=5e-08) 
show(p)


# In[ ]:


#make a QQ plot for TC anc0
p = hl.plot.qq(mt.results.TC.p_value[2], title="QQ plot, TC, anc0")
show(p)


# In[ ]:


#make a QQ plot for TC anc1
p = hl.plot.qq(mt.results.TC.p_value[3], title="QQ plot, TC, anc1")
show(p)

