![](images/tractor_icon.png)

# TRACTOR - Local Ancestry Aware GWAS

**Current version: 1.1.0**

Tractor is a specialized tool designed to enhance Genome-Wide Association Studies (GWAS) for diverse cohorts by addressing challenges associated with analyzing admixed populations. Admixed populations are often excluded from genomic studies due to concerns about how to properly account for their complex ancestry.

Tractor facilitates the inclusion of admixed individuals in association studies by leveraging local ancestry, allowing for finer resolution in controlling for ancestry in GWAS, and empowering identification of ancestry-specific loci.

## Classic GWAS vs. TRACTOR GWAS
Unlike traditional GWAS methods, Tractor requires local ancestry estimates in its analyses. It employs a multi-step approach involving phasing, local ancestry inference, and regression analysis with ancestral allele dosages. This method aims to improve the accuracy of association analyses in cohorts with diverse ancestries, overcoming issues such as population stratification and variable linkage disequilibrium patterns.

## Contents
* [Setup Conda environment](#setup-conda-environment)
* [Steps for Running Tractor Locally](#steps-for-running-tractor-locally)
  * [Optional Step: Recovering Haplotypes Disrupted by Statistical Phasing](#optional-step-recovering-haplotypes-disrupted-by-statistical-phasing)
  * [Extracting Tracts and Ancestry Dosages](#extracting-tracts-and-ancestry-dosages)
  * [Running Tractor](#running-tractor)
* [Steps for Running Tractor on Hail](#steps-for-running-tractor-on-hail)
* [Output Files](#output-files)
* [License](#license)
* [Cite this article](#cite-this-article)

## Setup Conda environment

We recommend creating a Conda environment to run Tractor locally.
```bash
conda env create -f conda_py3_tractor.yml
conda activate py3_tractor
```

[Contents](#contents)

## Steps for Running Tractor Locally

#### IMPORTANT: Ensure your genotype data is phased (VCF file) and local ancestry is inferred for the following steps. Refer to our [Tractor tutorial](https://atkinson-lab.github.io/Tractor-tutorial/) for initial setup instructions.

All scripts desribed in the following steps are available in the [`scripts`](https://github.com/Atkinson-Lab/Tractor-New/tree/main/scripts) directory, and Hail implementation is present in the [`ipynbs`](https://github.com/Atkinson-Lab/Tractor-New/tree/main/ipynbs) directory

### Optional Step: Recovering Haplotypes Disrupted by Statistical Phasing

Statistical phasing can lead to switch errors as described in [Fig. 1](https://www.nature.com/articles/s41588-020-00766-y/figures/1) of the Tractor publication.
For this purpose, we have written two scripts, `unkink_2way_mspfile.py` and `unkink_2way_genofile.py`. These scripts help recover disrupted tracts from the **MSP file and VCF file**, rectifying errors, and outputs an unkinked VCF file that can be used for subsequent steps. Currently they are implemented for two-way admixed popuations only.
- `unkink_2way_mspfile.py`
  ```
  --msp              Path stem to MSP file, not including ".msp.tsv". (Must end in .msp.tsv)
  ```
- **Output File:**
  - The output file \*.switches.txt includes information on windows from the MSP file that needs to be switched. 
  - This file will serve as an input to `unkink_2way_genofile.py`
- `unkink_2way_genofile.py`
  ```
  --switches         Path to *.switches.txt, which includes info on windows to be switched
  --genofile         Path stem to input VCF with phased genotypes, not including .vcf suffix
  ```
- **Output File:**
  - The output file \*.unkinked.vcf will include recovered haplotypes

[Contents](#contents)

### Extracting Tracts and Ancestry Dosages

Simultaneously extract risk allele and local ancestry information, a prerequisite for running Tractor GWAS. The scripts output risk allele by ancestry dosages and haplotype counts for the input VCF files. A file of each of these is generated for each ancestry component.
- **Note that the input VCF file must be the phased file on which local ancestry was called.**
- Running `ExtractTracts.py` requires the **input MSP and VCF file**, and the number of ancestral populations within the VCF file. This script outputs the dosage and hapcount files required for running Tractor using `run_tractor.R`.
  - `ExtractTracts.py`:
    ```
    --msp              Path stem to MSP file (do not include *.msp.tsv suffix)
    --vcf-prefix       Path stem to input VCF file with phased genotypes (do not include *.vcf or
                       *.vcf.gz suffix. If zipped, used --zipped argument)
    --num-ancs         Number of ancestral populations within the VCF/MSP file.
    --zipped           If input VCF file is gzipped
    --zip-output       If all output files need to be Gzip compressed. (not BGZF-compression, recommended for very large files)
    --output-path      Optional output path for files and file prefix, e.g. ~/test_data/test1
    ```

    **Example run for 3-way admixed dataset:**
    - If `dataset_qc_phased.msp.tsv` and `dataset_qc_phased.vcf` exists, we can run the following --
    ```bash
    python3 ExtractTracts.py \
    --msp dataset_qc_phased \
    --vcf-prefix dataset_qc_phased \
    --num-ancs 3 \
    --output-path ${output_dir}/output_prefix
    ```
    - Note: With some LAI algorithms, `*.msp` file is generated. In this case, please create a symbolic link to `*.msp.tsv` as the `ExtractTracts.py` expects MSP files to have `*.msp.tsv` suffix.

- If you used **[FLARE](https://github.com/browning-lab/flare) for LAI** and there are no MSP files, please await update as we are working towards a script that will help extract the necessary dosage and haplotype count files from FLARE VCF output.

- **Output Files:**
  - `ExtractTracts.py` will generate two files (`*.dosage.txt`, `*.hapcount.txt`) per ancestry.

[Contents](#contents)

### Running Tractor

- The Tractor code runs in R, and all required library packages should be installed within the Conda environment.
- Arguments:
  ```
  --hapdose         Prefix of hapcount and dosage files generated
  --phe             Phenotype file; 1st column sample ID, 2nd column phenotype,
                    other columns will be treated as covariates. Missing data is allowed.
  --method          "linear" or "logistic"
  --out             Output file name for ancestry-specific summary statistics
  ```

  **Example run:**
  ```
  ${script_path}/run_tractor.R \
  --hapdose dataset_qc_phased \
  --phe dataset_qc_pheno_covars.txt \
  --method logistic \
  --out dataset_qc_phased_sumstats
  ```

- **Output Files:**
  - Tractor is a local-ancestry aware GWAS that offers ancestry-specific summary statistics.
  - The number of columns would depend on the number of ancestries within the study. Here is a description of the columns of 2-way admixed dataset:
    ```
    CHROM:              Chromosome 
    POS:                Position 
    ID:                 SNP ID
    REF:                Reference allele 
    ALT:                Alternate allele 
    AF_anc0:            Allele frequency for anc0; sum(dosage)/sum(local ancestry)
    AF_anc1:            Allele frequency for anc1; sum(dosage)/sum(local ancestry)
    LAprop_anc0:        Local ancestry proportion for anc0; sum(local ancestry)/2 * sample size
    LAprop_anc1:        Local ancestry proportion for anc1; sum(local ancestry)/2 * sample size
    LAeff_anc0:         Effect size for the local ancestry term (X1 term in Tractor)
    LApval_anc0:        p value for the local ancestry term (X1 term in Tractor)
    Geff_anc0:          Effect size for alternate alleles that are interited from anc0
    Geff_anc1:          Effect size for alternate alleles that are interited from anc1
    Gpval_anc0:         p value for alternate alleles that are interited from anc0
    Gpval_anc1:         p value for alternate alleles that are interited from anc1
    ```
  - Gpval (Genotype p-value) columns can be used for generating ancestry-specific Manhattan plots.

[Contents](#contents)

## Steps for Running Tractor on Hail
- Hail implementation of the pipeline is described in [`hail_example_tractor_gwas.ipynb`](https://github.com/Atkinson-Lab/Tractor-New/blob/main/ipynbs/hail_example_tractor_gwas.ipynb).

[Contents](#contents)

## License
The Tractor program is licensed under the MIT License. You may obtain a copy of the License [here](https://github.com/Atkinson-Lab/Tractor-New/blob/main/LICENSE).

[Contents](#contents)

## Cite this article

The methodology and utility of Tractor are more fully described in our manuscript. If you use Tractor in your research, please cite the following article:

> Atkinson, E.G., Maihofer, A.X., Kanai, M. et al. Tractor uses local ancestry to enable the inclusion of admixed individuals in GWAS and to boost power. Nat Genet 53, 195–204 (2021). [Link](https://doi.org/10.1038/s41588-020-00766-y)

For any inquiries, you can contact Elizabeth G. Atkinson at [elizabeth.atkinson@bcm.edu](mailto:elizabeth.atkinson@bcm.edu).

[Contents](#contents)


