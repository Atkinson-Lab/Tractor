# The Tractor pipeline

A pipeline for improving the long-range phasing of genotype files informed by local ancestry deconvolution, and extracting the separate ancestral segments from admixed individuals. The pipeline is currently designed assuming RFmix_v2 formatted output and phased VCF format genotype input.

The pipeline consists of 3 scripts, the first two being optional:
1) Detect and unkink strand flips using RFmix .msp.tsv output. A strand flip is defined as a switch within a 1 cM window to opposite strands of at least one component ancestries (currently implemented in a 2-way admixed setting) conditioned on heterozygous ancestral dosage.
This step is implememted with the script *UnkinkMSPfile.py*, which tracks switch locations and unkinks your ancestry call RFmix output file. This recovers long range tracts which are disrupted by statistical phasing.

Example usage: 
```python UnkinkMSPfile.py --msp FILENAME_STEM```

2) Unkink the phased genotype file (VCF format) that was fed into RFmix. This improves long-range phasing accuracy in a manner informed by the local ancestry tract calls. This step is implememted with the script *UnkinkGenofile.py* and expects as input the phased VCF file that was fed into RFmix and the 'switches' file generated with step 1, which is used to determine the positions that need to be flipped in the VCF file. 

Example usage: 
```python UnkinkGenofile.py --switches SWITCHES_FILE --genofile INPUT_VCF```

3) Extract the tracts that relate to each component ancestry into their own VCFs. These ancestry-deconvoluted files can then be analyzed alongside larger cohorts of the ancestral composition they pertain to. For example, an admixed African American individual could be run through the pipeline and separated into "European" and "African" chromosomal fragments to be considered alongside European and African cohort collections rather than being excluded from study due to population structure concerns.
This step is implemented with the script *ExtractTracts.py* and expects as input the VCF and RFmix msp files. You can input the unkinked files from the previous step if you chose to unkink. The script *ExtractTracts-3way.py* handles three-way admixed population scenarios.

Example usage: 
```python ExtractTracts.py --msp MSPFILE --vcf VCF```


Please cite our publication:
[insert pub]

Other usage notes:
A thorough description of Tractor and its potential uses can be found in our paper, listed above. We additionally have developed a GWAS application, implemented in joint analysis format in Hail. This is distributed in the form of a Jupyter notebook, the file XX.

Important: The success of Tractor relies on good local ancestry calls. Please ensure your LAI performance is highly accurate. Additionally, at the ends of chromosomes it is hard to infer local ancestry accurately, so you may wish to trim the ends or treat chromosome ends with caution.
