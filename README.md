## The Tractor pipeline

Thanks for your interest in Tractor!

Please check out the Wiki page of this repo (https://github.com/eatkinson/Tractor/wiki) for more details about each step and descriptions of implementation.

The method and utility of Tractor is more fully described in our preprint, "Tractor: A framework allowing for improved inclusion of admixed individuals in large-scale association studies", which is on bioRxiv here: https://www.biorxiv.org/content/10.1101/2020.05.17.100727v1.


In this repo, we provide scripts comprising the full Tractor pipeline for (1) recovering long-range haplotypes as informed by local ancestry deconvolution, (2) extracting ancestral segments and calculating ancestral minor allele dosages, and (3) running our local ancestry aware GWAS model.

The pipeline is currently designed assuming RFmix_v2 and phased VCF format genotype input and consists of separate scripts for maximum flexibility of use. These scripts can be run in sequence to do all the steps or, if long-range tracts will not impact the goals of the analysis, the user can skip right to ancestry dosage calculation and running our GWAS model.


**Important**: The success of Tractor relies on good local ancestry calls. Please ensure your LAI performance is highly accurate before running a Tractor GWAS. We recommend admixture simulations to test accuracy for your demographic use case, which can be done with a program such as https://github.com/williamslab/admix-simu. 
