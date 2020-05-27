## The Tractor pipeline

Thanks for your interest in Tractor!

The method and implementation is more fully described in our preprint, "Tractor: A framework allowing for improved inclusion of admixed individuals in large-scale association studies", which is on bioRxiv here.

In this repo, we provide scripts comprising the full Tractor pipeline for (1) recovering long-range haplotypes as informed by local ancestry deconvolution, (2) extracting ancestral segments from admixed individuals, (3) calculating ancestral minor allele dosages and (4) running our local ancestry aware GWAS model. The pipeline is currently designed assuming RFmix_v2 and phased VCF format genotype input.

The pipeline consists of separate scripts for maximum flexibility of utility. These scripts can be run in sequence to do all the steps or, if long-range tracts will not impact the goals of the analysis, the user can skip right to ancestry dosage calculation and running our GWAS model.

Please check out the Wiki page of this repo for details about each step.
