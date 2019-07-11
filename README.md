# The Tractor pipeline

A pipeline for improving the long-range phasing of genotype files informed by local ancestry deconvolution, and extracting the separate ancestral segments from admixed individuals. The pipeline is currently designed assuming RFmix_v2 formatted output and VCF format genotype input.

The pipeline consists of scripts for each of 3 steps:
1) Detect and unkink strand flips using RFmix .msp.tsv output. A strand flip is defined as a switch within a 1 cM window to opposite strands of at least one component ancestries (currently implemented in a 2-way admixed setting) conditioned on heterozygous ancestral dosage.
This step is implememted with the script *UnkinkMSPfile-Flags-switchpoints.py*, which unkinks your ancestry call RFmix output file and logs the location of switch points. 

Example usage: 
```python UnkinkMSPfile-Flags-switchpoints.py --msp FILENAME_STEM```

2) Unkink the phased genotype file (VCF format) that was fed into RFmix. This improves long-range phasing accuracy in a manner informed by the local ancestry tract calls. This step is implememted with the script *UnkinkGenofile-Flags.py* and expects as input the phased VCF file that was fed into RFmix and the 'switches' file generated with step 1, which is used to determine the positions that need to be flipped in the VCF file. 

Example usage: 
```python UnkinkGenofile-Flags.py --switches SWITCHES_FILE --genofile INPUT_VCF```

3) Extract the tracts that relate to each component ancestry into their own VCFs. These ancestry-deconvoluted files can then be analyzed alongside larger cohorts of the ancestral composition they pertain to. For example, an admixed African American individual could be run through the pipeline and separated into "European" and "African" chromosomal fragments to be considered alongside European and African cohort collections rather than being excluded from study due to population structure concerns.
This step is implemented with the script *ExtractTracts-Flags.py* and expects as input the now-unkinked VCF and RFmix msp files.

Example usage: 
```python ExtractTracts-Flags.py --msp UNKINKED_MSPFILE --vcf UNKINKED_VCF```
