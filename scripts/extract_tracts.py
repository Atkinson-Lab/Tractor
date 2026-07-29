#!/usr/bin/env python3
# Tractor script to extract out ancestry segments from each reference population from admixed samples.
# This script has been modified to work with RFMix2 and G-Nomix output, which includes ancestry calls
# in the *.msp file and phased genotypic data in *.vcf (or *.vcf.gz) file.

# Original code by Dr. Elizabeth G. Atkinson (elizabeth.atkinson@bcm.edu)
__author__ = "elizabeth.atkinson@bcm.edu"
__version__ = "1.2.0"

# Changelog:
# Oct 14, 2023 (v1.1.0)
# Modified input arguments, --vcf and --msp take complete file names now.
# --zip-output is now --compress-output. Output of ancestry specific VCF files
# is now optional, with --output-vcf argument.
# Feb 20, 2024 (v1.2.0)
# Discrepancy of output between extract_tracts.py (v1.1.2) and ExtractTracts.py
# was due to an unintended "if output_vcf:" in code when computing dosage. This
# has been removed for v1.2.0.

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
Usage:
    extract_tracts.py --vcf <Required. Path to VCF file (*.vcf or *.vcf.gz)>
                      --msp <Required. Path to MSP file (*.msp or *.msp.tsv)>
                      --num-ancs <Required. Number of ancestral populations in the VCF file>
                      --output-dir <Optional. Path to output directory (default: same directory as VCF file)>
                      --output-vcf <Optional. Default: False. Output ancestry-specific VCF file>
                      --compress-output <Optional. Compress output dosage, hapcount, and VCF files>
Options:
    --vcf                Path to VCF file (*.vcf or *.vcf.gz)
    --msp                Path to MSP file (*.msp or *.msp.tsv)
    --num-ancs           Number of ancestral populations within the VCF file.
    --output-dir         Path to the output directory (default: same as VCF file). Must already exist.
    --output-vcf         Whether to output the ancestry-specific VCF file (default: False).
    --compress-output    Compress output dosage, hapcount, and VCF files (default: False).
"""


def parse_vcf_filepath(vcf_file_path):

    if not os.path.exists(vcf_file_path):
        raise FileNotFoundError(f"File not found: {vcf_file_path}")

    if os.path.getsize(vcf_file_path) == 0:
        err_msg = f"The provided vcf file, {vcf_file_path}, was empty"
        logger.critical(err_msg)
        raise ValueError(err_msg)

    zipped = False
    path, filename = os.path.split(vcf_file_path)

    # check if the file is compressed (ends with '.gz')
    if filename.endswith(".vcf.gz"):
        zipped = True
        prefix = os.path.splitext(os.path.splitext(filename)[0])[
            0
        ]  # remove '.vcf.gz' extension
    elif filename.endswith(".vcf"):
        prefix = os.path.splitext(filename)[0]  # remove '.vcf' extension
    else:
        # File has an unexpected extension
        raise ValueError(
            f"Unexpected file extension: {filename}. VCF file must end in .vcf or .vcf.gz"
        )

    return {"path": path, "prefix": prefix, "zipped": zipped}


def extract_tracts(
    vcf=str,
    msp=str,
    num_ancs=int,
    output_dir=None,
    output_vcf=None,
    compress_output=None,
):

    # Housekeeping checks.
    if not os.path.exists(vcf):
        raise ValueError(f"The path '{vcf}' does not exist.")
    if not os.path.exists(msp):
        raise ValueError(f"The path '{msp}' does not exist.")
    if not isinstance(num_ancs, int):
        raise ValueError("The 'num_ancs' must be an integer.")
    if output_dir is not None:
        if not os.path.exists(output_dir):
            raise ValueError(f"The output directory '{output_dir}' does not exist.")
    if not isinstance(output_vcf, bool):
        raise ValueError("The 'output_vcf' must be a boolean value (True or False).")
    if not isinstance(compress_output, bool):
        raise ValueError(
            "The 'compress_output' must be a boolean value (True or False)."
        )
        if compress_output == True:
            logger.info("Files will be gzip compressed (not bgzip)")

    vcf_info = parse_vcf_filepath(vcf)
    vcf_path, vcf_prefix, zipped = (
        vcf_info["path"],
        vcf_info["prefix"],
        vcf_info["zipped"],
    )

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
        output_files[f"dos{i}"] = f"{output_path}.anc{i}.dosage.txt{file_extension}"
        output_files[f"ancdos{i}"] = (
            f"{output_path}.anc{i}.hapcount.txt{file_extension}"
        )
        if output_vcf:
            output_files[f"vcf{i}"] = f"{output_path}.anc{i}.vcf{file_extension}"

    with open(msp) as mspfile, (
        gzip.open(vcf, "rt") if zipped else open(vcf)
    ) as vcf, contextlib.ExitStack() as stack:
        files = {
            fname: stack.enter_context(
                gzip.open(output_file, "wt")
                if compress_output
                else open(output_file, "w")
            )
            for fname, output_file in output_files.items()
        }
        vcf_header = ""
        window = ("", 0, 0)  # initialize the current window to check

        logger.info("Iterating through VCF file")
        for line in vcf:
            skip_line = False
            if line.startswith("##"):
                if output_vcf:
                    vcf_header += line
                continue
            elif line.startswith("#"):
                if output_vcf:
                    vcf_header += line
                anc_header = "\t".join(
                    [
                        line.strip("# ").split("\t", 9)[item]
                        for item in [0, 1, 2, 3, 4, 9]
                    ]
                )

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
                pos = int(row[1])
                output_lines = {}
                pop_genos = {}

                for i in range(num_ancs):  # Write entries into each output files' list
                    output_lines[f"dos{i}"] = dos_anc_out
                    output_lines[f"ancdos{i}"] = dos_anc_out
                    if output_vcf:
                        output_lines[f"vcf{i}"] = vcf_out

                # optimized for quicker runtime - only move to next line when out of the current msp window
                # saves the current line until out of the window, then checks next line. VCF and window switches file should be in incremental order.

                while not (row[0] == window[0] and (window[1] <= pos < window[2])):
                    if row[0] == window[0] and window[1] > pos:
                        skip_line = True  # Skip VCF line
                        break  # Break out of msp scanning because the VCF position is still before the windows start
                    ancs = mspfile.readline()
                    if ancs.startswith("#"):  # skip the header lines
                        continue
                    if not ancs:
                        break  # when get to the end of the msp file, stop
                    # chm, spos, epos, sgpos, egpos, nsnps, calls
                    ancs_entry = ancs.strip().split("\t", 6)
                    calls = ancs_entry[6].split("\t")
                    window = (ancs_entry[0], int(ancs_entry[1]), int(ancs_entry[2]))
                    if row[0] == window[0] and window[1] > pos:
                        skip_line = True  # Skip VCF line
                        break  # Break out of msp scanning because the VCF position is before the windows start
                if skip_line:
                    logger.info(
                        f"VCF position, {pos} is not in an msp window, skipping site"
                    )
                    continue  # Skip VCF site and continue to next site

                # index by the number of individuals in the VCF file, should be half the number in the calls file
                for i, geno in enumerate(genos):
                    geno_parts = geno.split(":")[0].split(
                        "|"
                    )  # assert incase eagle leaves some genos unphased
                    geno_a, geno_b = map(str, geno_parts[:2])
                    call_a = str(calls[2 * i])
                    call_b = str(calls[2 * i + 1])

                    counts = {anc: 0 for anc in range(num_ancs)}
                    anc_counts = {anc: 0 for anc in range(num_ancs)}
                    for j in range(num_ancs):
                        if output_vcf:
                            pop_genos[j] = ""
                        if call_a == str(j):
                            if output_vcf:
                                pop_genos[j] += geno_a
                            anc_counts[j] += 1
                            if geno_a == "1":
                                counts[j] += 1
                        else:
                            if output_vcf:
                                pop_genos[j] += "."

                        if call_b == str(j):
                            if output_vcf:
                                pop_genos[j] += "|" + geno_b
                            anc_counts[j] += 1
                            if geno_b == "1":
                                counts[j] += 1
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

    logger.info("Finished extracting tracts per %d ancestries", num_ancs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vcf",
        required=True,
        help="Path to phased genotypes, VCF file (*.vcf or *.vcf.gz)",
    )
    parser.add_argument(
        "--msp",
        required=True,
        help="Path to ancestry calls, MSP file (*.msp or *.msp.tsv)",
    )
    parser.add_argument(
        "--num-ancs",
        type=int,
        required=True,
        help="Number of ancestral populations within the VCF file.",
    )
    parser.add_argument(
        "--output-dir", help="Directory for output files. Directory must already exist."
    )
    parser.add_argument(
        "--output-vcf",
        help="Boolean. If ancestry-specific VCF files need to be generated",
        action="store_true",
    )
    parser.add_argument(
        "--compress-output",
        help="gzip (not bgzip) all output files",
        action="store_true",
    )

    args = parser.parse_args()
    extract_tracts(**vars(args))
