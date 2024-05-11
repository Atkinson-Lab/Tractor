#!/usr/bin/env python3

# Tractor script to extract out ancestry segments from each reference population from admixed samples.
# The original ExtractTracts.py has been modified to work with FLARE output wherein ancestry calls
# are saved within the VCF file itself, and no MSP file is generated.
# Several changes were made, including addition of several arguments, but the basic algorithm for
# creating dosage and hapcount files remain.

# Original code by Dr. Elizabeth G. Atkinson (elizabeth.atkinson@bcm.edu)
__author__ = "elizabeth.atkinson@bcm.edu"
__version__ = "1.0.1"

# Modified by Nirav N. Shah (nirav.shah@bcm.edu)
# Changelog:
# Oct 14, 2023 (v1.0.0)
# Compatible with FLARE v0.3.0 output.

import argparse
import contextlib
import gzip
import logging
import re
import os

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

USAGE = """
This script generates ancestry-specific dosage and hapcount files from an input VCF,
similar to or identical to FLARE's output VCF file, which includes genotype and
ancestry calls (tested on FLARE v0.3.0 output).

Usage:
    extract_tracts_flare.py --vcf <Required. Path to VCF file (*.vcf or *.vcf.gz)>
                            --num-ancs <Required. Number of ancestral populations in the VCF file>
                            --output-dir <Optional. Path to output directory (default: same directory as VCF file)>
                            --output-vcf <Optional. Default: False. Output ancestry-specific VCF file>
                            --compress-output <Optional. Compress output dosage, hapcount, and VCF files>

Options:
    --vcf                   Path to VCF file (*.vcf or *.vcf.gz)
    --num-ancs              Number of ancestral populations within the VCF file.
    --output-dir            Path to the output directory (default: same as VCF file). Must already exist.
    --output-vcf            Whether to output the ancestry-specific VCF file (default: False).
    --compress-output       Compress output dosage, hapcount, and VCF files (default: False).

"""

def parse_vcf_filepath(vcf_file_path):

    if not os.path.exists(vcf_file_path):
        raise FileNotFoundError(f"File not found: {vcf_file_path}")

    zipped = False
    path, filename = os.path.split(vcf_file_path)

    # check if the file is compressed (ends with '.gz')
    if filename.endswith('.vcf.gz'):
        zipped = True
        prefix = os.path.splitext(os.path.splitext(filename)[0])[0] # remove '.vcf.gz' extension
    elif filename.endswith('.vcf'):
        prefix = os.path.splitext(filename)[0] # remove '.vcf' extension
    else:
        # File has an unexpected extension
        raise ValueError(f"Unexpected file extension: {filename}. VCF file must end in .vcf or .vcf.gz")

    return {
        'path': path,
        'prefix': prefix,
        'zipped': zipped
    }

def extract_tracts_flare(vcf=str, num_ancs=int, output_dir=None, output_vcf=None, compress_output=None):

    # Print version of the script
    logger.info(f"# This is version {__version__} of the extract_tracts_flare.py script.")
    
    # Housekeeping checks.
    if not os.path.exists(vcf):
        raise ValueError(f"The path '{vcf}' does not exist.")
    if not isinstance(num_ancs, int):
        raise ValueError("The 'num_ancs' must be an integer.")
    if output_dir is not None:
        if not os.path.exists(output_dir):
            raise ValueError(f"The output directory '{output_dir}' does not exist.")
    if not isinstance(output_vcf, bool):
        raise ValueError("The 'output_vcf' must be a boolean value (True or False).")
    if not isinstance(compress_output, bool):
        raise ValueError("The 'compress_output' must be a boolean value (True or False).")
        if compress_output==True:
            logger.info("Files will be gzip compressed (not bgzip)")
    
    vcf_info = parse_vcf_filepath(vcf)
    vcf_path, vcf_prefix, zipped = vcf_info['path'], vcf_info['prefix'], vcf_info['zipped']

    # Log version number
    logger.info("# Running script version      : %s", __version__)
    logger.info("# VCF File                    : %s", vcf)
    logger.info("# Prefix of output file names : %s", vcf_prefix)
    logger.info("# VCF File is compressed?     : %s", zipped)
    logger.info("# Number of Ancestries in VCF : %d", num_ancs)
    logger.info("# Output Directory            : %s", output_dir)
    
    output_files = {}
    output_path = f"{os.path.join(output_dir, vcf_prefix) if output_dir else os.path.join(vcf_path, vcf_prefix)}"

    file_extension = f"{'.gz' if compress_output else ''}"

    # Output files: vcf, dosage, haplotype count per passed ancestry
    logger.info("Creating output files for %d ancestries", num_ancs)
    for i in range(num_ancs):
        if output_vcf:
            output_files[f"vcf{i}"] = f"{output_path}.anc{i}.vcf{file_extension}" 
        output_files[f"dos{i}"] = f"{output_path}.anc{i}.dosage.txt{file_extension}"
        output_files[f"ancdos{i}"] = f"{output_path}.anc{i}.hapcount.txt{file_extension}"

    logger.info("Opening input and output files for reading and writing")
    with gzip.open(vcf, "rt") if zipped else open(vcf) as vcf, contextlib.ExitStack() as stack:
        files = {
            fname: stack.enter_context(gzip.open(output_file, "wt") if compress_output else open(output_file, "w"))
            for fname, output_file in output_files.items()
        }
        vcf_header = ""
        
        logger.info("Iterating through VCF file")
        for line in vcf:
            if line.startswith("##"):
                if output_vcf:
                    vcf_header += line 
                continue
            elif line.startswith("#"):
                if output_vcf:
                    vcf_header += line
                anc_header = "\t".join([ line.strip("# ").split("\t", 9)[item] for item in [0, 1, 2, 3, 4, 9] ])

                for i in range(num_ancs):  # Write header into each output vcf
                    files[f"dos{i}"].write(anc_header)
                    files[f"ancdos{i}"].write(anc_header)
                    if output_vcf:
                        files[f"vcf{i}"].write(vcf_header)
            
            if not line.startswith("#"):
                # Entry format is ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'genotypes']
                row = line.strip().split("\t", 9)

                # Grab fields needed for each file output
                row[8] = "GT"  # Update FORMAT field to only keep GT
                vcf_out = "\t".join(row[:9])
                dos_anc_out = "\t".join(row[:5])

                # split each strand geno call apart from each other
                genos = row[9].split("\t")
                output_lines = {}
                pop_genos = {}

                for i in range(num_ancs):  # Write entries into each output files' list
                    output_lines[f"dos{i}"] = dos_anc_out
                    output_lines[f"ancdos{i}"] = dos_anc_out
                    if output_vcf:
                        output_lines[f"vcf{i}"] = vcf_out
             
                # Genotype data is of format GT1|GT2:ANC1:ANC2, where GT1:GT2 is the genotype, 
                # and ANC1:ANC2 are the ancestry calls
                for i,geno in enumerate(genos):
                    geno_parts = re.split('[|:]',geno) 

                    geno_a, geno_b, call_a, call_b = map(str, geno_parts[:4])

                    counts     = {anc: 0 for anc in range(num_ancs)}
                    anc_counts = {anc: 0 for anc in range(num_ancs)}

                    for j in range(num_ancs):
                        if output_vcf:
                            pop_genos[j] = ""
                        if call_a == str(j):
                            if output_vcf:
                                pop_genos[j]+= geno_a
                            anc_counts[j]+=1
                            if geno_a == "1":
                                counts[j]+=1
                        else:
                            if output_vcf:
                                pop_genos[j] += "."
                        
                        if call_b == str(j):
                            if output_vcf:
                                pop_genos[j] += "|" + geno_b
                            anc_counts[j]+=1
                            if geno_b == "1":
                                counts[j]+=1
                        else:
                            if output_vcf:
                                pop_genos[j] += "|."
                        
                        output_lines[f"dos{j}"] += "\t" + str(counts[j])
                        output_lines[f"ancdos{j}"] += "\t" + str(anc_counts[j])
                        if output_vcf:
                            output_lines[f"vcf{j}"] += "\t" + pop_genos[j]

                for j in range(num_ancs):
                    output_lines[f"dos{j}"] += "\n"
                    output_lines[f"ancdos{j}"] += "\n"
                    files[f"dos{j}"].write(output_lines[f"dos{j}"])
                    files[f"ancdos{j}"].write(output_lines[f"ancdos{j}"])
                    if output_vcf:
                        output_lines[f"vcf{j}"] += "\n"
                        files[f"vcf{j}"].write(output_lines[f"vcf{j}"])
    
    logger.info("Finished extracting tracts for %d ancestries.", num_ancs)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf",required=True,help="Path to VCF file (*.vcf or *.vcf.gz)")
    parser.add_argument("--num-ancs",type=int,required=True,
                        help="Number of ancestral populations within the VCF file.")
    parser.add_argument("--output-dir", help="Directory for output files. Directory must already exist.")
    parser.add_argument("--output-vcf", help="Boolean. If ancestry-specific VCF files need to be generated",
                        action="store_true")
    parser.add_argument("--compress-output", help="gzip (not bgzip) all output files",
                        action="store_true")

    args = parser.parse_args()
    extract_tracts_flare(**vars(args))
    
    logger.info("# Run complete")

