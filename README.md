## The Tractor pipeline

Thanks for your interest in Tractor!

Please check out the Wiki page of this repo for more details about each step and descriptions of implementation.

The method and implementation is more fully described in our preprint, "Tractor: A framework allowing for improved inclusion of admixed individuals in large-scale association studies", which is on bioRxiv here: https://www.biorxiv.org/content/10.1101/2020.05.17.100727v1.

In this repo, we provide scripts comprising the full Tractor pipeline for (1) recovering long-range haplotypes as informed by local ancestry deconvolution, (2) extracting ancestral segments from admixed individuals, (3) calculating ancestral minor allele dosages and (4) running our local ancestry aware GWAS model. The pipeline is currently designed assuming RFmix_v2 and phased VCF format genotype input.

The Tractor pipeline consists of separate scripts for maximum flexibility of utility. These scripts can be run in sequence to do all the steps or, if recovering long-range tracts will not impact the goals of the analysis, the user can skip right to ancestry dosage calculation and conducting an ancestry-aware association study.


**Important**: The success of Tractor relies on good local ancestry calls. Please ensure your LAI performance is highly accurate before running a Tractor GWAS. We recommend admixture simulations to test accuracy for your demographic use case, which can be implemented with a program such as https://github.com/williamslab/admix-simu. Additionally, at the ends of chromosomes it is hard to infer local ancestry accurately, so you may wish to trim the ends or treat chromosome ends with caution.
