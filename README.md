![](images/tractor_icon.png)

## NEW!!! Current Version: v1.4.0 (released May 10, 2024)
- Added support for compressed (gz) hapcount/dosage and phenotype files.
- Improved file reading efficiency by implementing fread in chunks, mitigating memory errors.
- Implemented parallel processing for regression, resulting in significant speed improvements with multi-core systems.
- Enhanced flexibility in organizing phenotype files:
   - Users can specify sample ID column (`--sampleidcol`), phenotype ID column (`--phenocol`), and covariate column list (`--covarcollist`)
- Updated output summary statistics to include SE and t-val, with column names adjusted to adhere to GWAS standards.

# TRACTOR - Local Ancestry Aware GWAS

Tractor is a specialized tool designed to enhance Genome-Wide Association Studies (GWAS) for diverse cohorts by addressing challenges associated with analyzing admixed populations. Admixed populations are often excluded from genomic studies due to concerns about how to properly account for their complex ancestry.

Tractor facilitates the inclusion of admixed individuals in association studies by leveraging local ancestry, allowing for finer resolution in controlling for ancestry in GWAS, and empowering identification of ancestry-specific loci.

## Classic GWAS vs. TRACTOR GWAS
Unlike traditional GWAS methods, Tractor requires local ancestry estimates in its analyses. It employs a multi-step approach involving phasing, local ancestry inference, and regression analysis with ancestral allele dosages. This method aims to improve the accuracy of association analyses in cohorts with diverse ancestries, overcoming issues such as population stratification and variable linkage disequilibrium patterns.



## Contents
* [Setup Conda environment](#setup-conda-environment)
* [Steps for Running Tractor Locally](#steps-for-running-tractor-locally)
  * [Step 0 \[Optional\]: Recovering Haplotypes Disrupted by Statistical Phasing](#step-0-optional-recovering-haplotypes-disrupted-by-statistical-phasing)
  * [Step 1: Extracting Tracts and Ancestry Dosages](#step-1-extracting-tracts-and-ancestry-dosages)
  * [Step 2: Running Tractor](#step-2-running-tractor) **(Tractor v1.4.0 released with additional functionalities)**
* [Output Files (Running Tractor)](#output-files-running-tractor)
* [Steps for Running Tractor on Hail](#steps-for-running-tractor-on-hail)
* [License](#license)
* [Cite this article](#cite-this-article)

## Setup Conda environment

We recommend creating a Conda environment to run Tractor locally. This will install the necessary Python 3 and R dependencies required by the scripts.
```bash
conda env create -f conda_py3_tractor.yml
conda activate py3_tractor
```

[Contents](#contents)

## Steps for Running Tractor Locally

#### IMPORTANT: Ensure your genotype data is phased (VCF file) and local ancestry is inferred for the following steps. Refer to our [Tractor tutorial](https://atkinson-lab.github.io/Tractor-tutorial/) for initial setup instructions.

All scripts desribed in the following steps are available in the [`scripts`](https://github.com/Atkinson-Lab/Tractor-New/tree/main/scripts) directory, and Hail implementation is present in the [`ipynbs`](https://github.com/Atkinson-Lab/Tractor-New/tree/main/ipynbs) directory

### Step 0 [Optional]: Recovering Haplotypes Disrupted by Statistical Phasing

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

### Step 1: Extracting Tracts and Ancestry Dosages

Simultaneously extract risk allele and local ancestry information, a prerequisite for running Tractor GWAS. The scripts output risk allele by ancestry dosages and haplotype counts for the input VCF files. A file of each of these is generated for each ancestry component.
- **Note that the input VCF file must be the phased file on which local ancestry was called.**
- Running `extract_tracts.py` requires the **input MSP and VCF file**, and the number of ancestral populations within the VCF file. This script outputs the dosage and hapcount files required for running Tractor.
  - `extract_tracts.py`:
    ```
    --vcf              Path to VCF file (*.vcf or *.vcf.gz)
    --msp              Path to MSP file (*.msp or *.msp.tsv)
    --num-ancs         Number of ancestral populations within the VCF file.
    --output-dir       Path to the output directory (default: same as VCF file).
                       Directory must already exist.
    --output-vcf       Whether to output the ancestry-specific VCF file (default: False).
    --compress-output  Compress output dosage, hapcount, and VCF files (default: False).
    ```

    **Example run for 3-way admixed dataset:**
    ```bash
    python3 extract_tracts.py \
    --vcf dataset_qc_phased.vcf \
    --msp dataset_qc_phased.msp \
    --num-ancs 3 \
    --output-dir ${output_dir}
    ```

- If you used **[FLARE](https://github.com/browning-lab/flare) for LAI** and there are no MSP files, use the `extract_tracts_flare.py` script instead.
  - `extract_tracts_flare.py`:
    ```
    --vcf              Path to VCF file (*.vcf or *.vcf.gz)
    --num-ancs         Number of ancestral populations within the VCF file.
    --output-dir       Path to the output directory (default: same as VCF file).
                       Directory must already exist.
    --output-vcf       Whether to output the ancestry-specific VCF file (default: False).
    --compress-output  Compress output dosage, hapcount, and VCF files (default: False).
    ```

    **Example run for 3-way admixed dataset (with FLARE LAI):**
    ```bash
    python3 extract_tracts_flare.py \
    --vcf dataset_qc_phased.vcf \
    --num-ancs 3 \
    --output-dir ${output_dir}
    ```

- **Output Files:**
  - Both scripts, `extract_tracts.py` and `extract_tracts_flare.py` will generate two files (`*.dosage.txt`, `*.hapcount.txt`) per ancestry.
  - Additional ancestry-specific VCF files may be generated if `--output-vcf` argument is provided. This file is not required for running Tractor, but might be needed for painting ancestry-specific individual karyograms for personal research purposes.
  - If files are compressed using the option `--compress-output`, they will undergo standard GZ compression, not BGZF compression (as expected by bcftools)

[Contents](#contents)

### Step 2: Running Tractor

The Tractor code runs in R, and to make sure the script works, you'll need to install the following libraries. Your conda environment should handle these installations by default.
```
install.packages('optparse')
install.packages('data.table')
install.packages('R.utils')
install.packages('dplyr')
install.packages('doParallel')
```

**Arguments:**
```
--hapdose       [Mandatory] Prefix for hapcount and dosage files.
                    E.g. If you have the following files:
                         filename.anc0.dosage.txt filename.anc0.hapcount.txt
                         filename.anc1.dosage.txt filename.anc1.hapcount.txt
                    use "--hapdose filename".
--phenofile     [Mandatory] Path to the file containing phenotype and covariate data. 
                    Default assumptions: Sample ID column: "IID" or "#IID", Phenotype column: "y".
                    If different column names are used, refer to --sampleidcol and --phenocol arguments.
                    All covariates MUST be included using --covarcollist.
--covarcollist  [Mandatory] Specify column names of covariates in the --phenofile.
                    Only listed columns will be included as covariates.
                    Separate multiple covariates with commas.
                    E.g. --covarcollist age,sex,PC1,PC2.
                    To exclude covariates, specify "--covarcollist none".
--method        [Mandatory] Specify the method to be used: <linear> or <logistic>.
--output        [Mandatory] File name for summary statistics output.
                    E.g. /path/to/file/output_sumstats.txt

--sampleidcol   [Optional] Specify sample ID column name in the --phenofile.
                    Default: "IID" or "#IID"
--phenocol      [Optional] Specify phenotype column name in the --phenofile.
                    Default: "y"
--chunksize     [Optional] Number of rows to read at once from hapcount and dosage files.
                    Use smaller values for lower memory usage.
                    Note: Higher chunksize speeds up streaming but requires more memory.
                    If out-of-memory errors occur, try increasing memory or
                    reducing --chunksize or --nthreads.
                    Default: 10000
--nthreads      [Optional] Specify number of threads to use.
                    Increasing threads can speed up processing but may increase memory usage.
                    Default: 1
--totallines    [Optional] Specify total number of lines in hapcount/dosage files (wc -l *.hapcount.txt).
                    If not provided, it will be calculated internally (recommended).
                    Exercise caution: if --totallines is smaller than the actual lines in the files, 
                    only a subset of data will be analyzed. If larger than the actual lines in the files,
                    an error will occur. Both scenarios are discouraged.
```

**Example Run (with Mandatory Arguments)**

- The latest Tractor v1.4.0 update introduces changes to default arguments, enhancing versatility and applicability.
- To run Tractor with the default assumptions, only 5 arguments are required.
- Ensure all covariates are specified using the `--covarcollist` flag.
```
run_tractor.R \
--hapdose /path/to/file/tmp1 \
--phenofile /path/to/file/dataset_qc_pheno_covars.txt \
--covarcollist age,sex,PC1,PC2,PC3,PC4,PC5 \
--method linear \
--output /path/to/results/test1.txt
```

**Example Run (with Optional Arguments)**

- In real-world scenarios, datasets may vary in size and default assumptions may not apply. Tractor accommodates these scenarios with optional arguments.
- Assuming a phenotype file with columns: PC1, PC2, PC3, PC4, PC5, age, sex, pheno1, pheno2, pheno3, sample_id, users can perform GWAS across different phenotypes.
- Below is an example to run Tractor GWAS for the **pheno1** phenotype:
```
run_tractor.R \
--hapdose /path/to/file/tmp1 \
--phenofile /path/to/file/dataset_qc_pheno_covars.txt \
--covarcollist age,sex,PC1,PC2,PC3,PC4,PC5 \
--method linear \
--output /path/to/results/test1.txt \
--sampleidcol sample_id \
--phenocol pheno1
```

- Users can utilize multi-threading for improved performance, and control file reading with chunking to avoid memory errors with extremely large files.
- Ensure a balance between chunk size and thread count to optimize performance without encountering memory issues.
- Below is an example to run Tractor GWAS for the **pheno1** phenotype w/ multithreading (4 cpu) and larger chunksize (15000):
```
run_tractor.R \
--hapdose /path/to/file/tmp1 \
--phenofile /path/to/file/dataset_qc_pheno_covars.txt \
--covarcollist age,sex,PC1,PC2,PC3,PC4,PC5 \
--method linear \
--output /path/to/results/test1.txt \
--sampleidcol sample_id \
--phenocol pheno1 \
--chunksize 15000 \
--nthreads 4
```

## Output Files (Running Tractor)

Tractor generates ancestry-specific summary statistics, producing output files with varying column numbers based on the input number of ancestries.

All summary statistic files include:
* **Variant Information:**
  * CHR: Chromosome
  * POS: Position
  * ID: SNP ID
  * REF: Reference allele
  * ALT: Alternate allele
* **Sample Size:**
  * N: Total number of samples going into the model (after exclude NAs).
    * Note this number can vary from the number of samples present in hapcount/dosage files, as there may be samples with NAs in the phenotype file which are eventually skipped.
* **Allele Frequency (AF), Local Ancestry Proportion (LAprop), Effect Size (beta), p-value (pval), and t-value (tval):**
  * For each ancestry term (anc), there are 'n' sets of columns. For instance, if there are n=2 ancestries, expect 2 sets of columns for each of these parameters.
* **Local Ancestry (LA) Related Columns:**
  * LApval: p-value for the local ancestry term (X1 term in Tractor)
  * LAeff: Effect size for the local ancestry term (X1 term in Tractor)
  * For 'n' ancestry terms (anc), expect 'n-1' sets of these columns. For example, if there are n=2 ancestries, expect 1 set of columns for each of these parameters.

**Example Output File Structure**
```
CHR             Chromosome 
POS             Position 
ID              SNP ID
REF             Reference allele
ALT             Alternate allele
N               Total sample size
AF_anc0         Allele frequency for anc0; sum(dosage)/sum(local ancestry)
LAprop_anc0     Local ancestry proportion for anc0; sum(local ancestry)/2 * sample size
beta_anc0       Effect size for alternate alleles inherited from anc0
se_anc0         Standard error for effect size (beta_anc0)
pval_anc0       p-value for alternate alleles inherited from anc0 (NOT -log10(pvalues))
tval_anc0       t-value for anc0
...
LApval_anc0     p-value for the local ancestry term (X1 term in Tractor)
LAeff_anc0      Effect size for the local ancestry term (X1 term in Tractor)
...

```

[Contents](#contents)

## Steps for Running Tractor on Hail
- Hail implementation of the pipeline is described in [`hail_example_tractor_gwas.ipynb`](https://github.com/Atkinson-Lab/Tractor-New/blob/main/ipynbs/hail_example_tractor_gwas.ipynb).

## License
The Tractor program is licensed under the MIT License. You may obtain a copy of the License [here](https://github.com/Atkinson-Lab/Tractor-New/blob/main/LICENSE).

## Cite this article
The methodology and utility of Tractor are more fully described in our manuscript. If you use Tractor in your research, please cite the following article:

> Atkinson, E.G., Maihofer, A.X., Kanai, M. et al. Tractor uses local ancestry to enable the inclusion of admixed individuals in GWAS and to boost power. Nat Genet 53, 195–204 (2021). [Link](https://doi.org/10.1038/s41588-020-00766-y)

For any inquiries, you can contact Elizabeth G. Atkinson at [elizabeth.atkinson@bcm.edu](mailto:elizabeth.atkinson@bcm.edu).

[Contents](#contents)

