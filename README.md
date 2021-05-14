## The Tractor pipeline

Thanks for your interest in Tractor!

Please check out the Wiki page of this repo (https://github.com/eatkinson/Tractor/wiki) for more details about each step and descriptions of implementation.

The methodology and utility of Tractor is more fully described in our manuscript, "Tractor uses local ancestry to enable the inclusion of admixed individuals in GWAS and to boost power", which is was published in Nature Genetics here: https://www.nature.com/articles/s41588-020-00766-y. We ask that you cite this publication for work utilizing the Tractor software.


In this repo, we provide scripts comprising the full Tractor pipeline for (1) recovering long-range haplotypes as informed by local ancestry deconvolution, (2) extracting ancestral segments and calculating ancestral minor allele dosages, and (3) running our local ancestry aware GWAS model.

The pipeline is currently designed assuming RFmix_v2 and phased VCF format genotype input and consists of separate scripts for maximum flexibility of use. These scripts can be run in sequence to do all the steps or, if long-range tracts will not impact the goals of the analysis, the user can skip right to ancestry dosage calculation and running our GWAS model. GWAS steps can scale to an arbitrary user-specified number of ancestries.


**Important**: The success of Tractor relies on good local ancestry calls. Please ensure your LAI performance is highly accurate before running a Tractor GWAS. We recommend admixture simulations to test accuracy for your demographic use case, which can be done with a program such as https://github.com/williamslab/admix-simu. 
