This directory contains genomic chunks that can be used for phasing by SHAPEIT5

Each file contains at least 4 columns:
> Column 1: Chunk index
> Column 2: Chromosome/contig ID
> Column 3: Chunk coordinates (including left and right buffers i.e. includes flanking regions)
> COlumn 4: Chunk coordinates excluding buffer regions

* Column 3 is what is used to phase and the buffer regions overlap across sequential contigs for `phase_ligate` to find out phasing concordance within the buffer regions using Switch Error Rate (SER) metric, which SHAPEIT5 team recommends to be above 80 or 90. Please see [SHAPEIT5 documentation](https://odelaneau.github.io/shapeit5/docs/documentation/ligate/)


The b38 chunks for 4cM and 20cM have been obtained from [SHAPEIT5's resources on its GitHub](https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38)

The b37 chunks were generated using TGP Phase 3 (GRCh37) dataset, which was QC'd for high-quality variants, SNPs, biallelic variants and MAF of 0.05% before generating chunks using GLIMPSE2 to create genomic chunks. See GLIMPSE2's [docuemntaton here](https://odelaneau.github.io/GLIMPSE/docs/documentation/chunk/).

The `chunks_fullchromosome` directory is created if the user wants to to run phasing for whole chromosome at once. It is genomic coordinate-agnostic, and should work for both genome coordinates.


Here’s an improved version of your README, organized for clarity and ease of understanding:

---

# Genomic Chunks for Phasing with SHAPEIT5

This directory contains genomic chunks that are prepared for use with the phasing tool SHAPEIT5.

Each file in this directory contains at least four columns:

1. **Column 1: Chunk Index**: A unique identifier for each chunk.
2. **Column 2: Chromosome/Contig ID**: The identifier for the chromosome or contig to which the chunk belongs.
3. **Column 3: Chunk Coordinates (with Buffers)**: The genomic coordinates of the chunk, including the left and right buffer regions (flanking regions).
4. **Column 4: Chunk Coordinates (without Buffers)**: The genomic coordinates of the chunk, excluding the buffer regions.

- **Column 3** (Chunk Coordinates with Buffers) is used for phasing, with overlapping buffer regions across contigs allowing `phase_ligate` to assess phasing concordance using the Switch Error Rate (SER). [SHAPEIT5's `ligate` documentation](https://odelaneau.github.io/shapeit5/docs/documentation/ligate/) recomments an SER above 80-90% for optimal results.

## Resources

- The **b38 chunks** for 4cM and 20cM were sourced directly from [SHAPEIT5's GitHub resources](https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38). File names were standardized across directories, and alternative versions with a different chromosome naming convention were created.
- The **b37 chunks** were generated from the TGP Phase 3 (GRCh37) dataset, which underwent quality control for high-quality variants, SNPs, biallelic variants, and a minor allele frequency (MAF) of 0.05%. These chunks were created using GLIMPSE2; see the [GLIMPSE2's `chunk` documentation](https://odelaneau.github.io/GLIMPSE/docs/documentation/chunk/) for more details.
- The `chunks_fullchromosome` directory contains files for users who wish to perform phasing on entire chromosomes. These chunks are compatible with both b37 and b38 genome coordinates.

All directories include files named `chunk_chr${chr}.txt` and `chunk_chr${chr}_rechr.txt`, differing only in chromosome nomenclature. Files with the format `chunk_chr${chr}.txt` use chromosome IDs like 1, 2, 3, ..., 22, while `chunk_chr${chr}_rechr.txt` uses the format `chr1`, `chr2`, `chr3`, ..., `chr22`.
