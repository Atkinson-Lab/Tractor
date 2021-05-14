# Tractor script to extract out ancestry segments from each reference population from admixed samples.
__author__ = "egatkinson"


import argparse
import gzip
import contextlib

USAGE = """
ExtractTracts.py --msp  <an ancestral calls file prefix produced by RFmix version 2, do not include file extension .msp.tsv>
                --vcf <VCF file prefix, do not include file extensions (.vcf, .vcf.gz) off>
                --zipped <Whether VCF is gzipped (stored as True so do not use unless VCF is gzipped))>
                --num-ancs <Number of ancestral populations within the msp file>
"""

# input is expected to be a VCF file suffixed with .vcf and an ancestral calls file produced by RFmix version 2, suffixed with .msp.tsv
# output will be returned in the same location as input VCF file with same naming convention


def extract_tracts(msp: str, vcf_prefix: str, zipped: bool, num_ancs: int = 2):
    """
    Extract ancestry segments from reference populations in admixed samples.

    Parse phased VCF into N ancestry-specific VCFs, dosage count files, and haplotype count files using a msp file sample ancestry look up.

    :param msp: Prefix to msp file.
    :param vcf_prefix: Prefix to VCF file.
    :param zipped: Whether of not the VCF file is gzipped.
    :param num_ancs: The number of ancestries within the msp file
    """
    mspfile = f"{msp}.msp.tsv"
    vcf = f"{vcf_prefix}.vcf.gz" if zipped else f"{vcf_prefix}.vcf"
    output_files = {}

    # Output files: vcf, dosage, haplotype count per passed ancestry
    for i in range(num_ancs):
        output_files[f"vcf{i}"] = f"{vcf_prefix}.anc{i}.vcf"
        output_files[f"dos{i}"] = f"{vcf_prefix}.anc{i}.dosage.txt"
        output_files[f"ancdos{i}"] = f"{vcf_prefix}.anc{i}.hapcount.txt"

    with open(mspfile) as mspfile, gzip.open(vcf, "rt") if zipped else open(
        vcf
    ) as vcf, contextlib.ExitStack() as stack:
        files = {
            fname: stack.enter_context(open(output_file, "w"))
            for fname, output_file in output_files.items()
        }
        vcf_header = ""
        window = ("", 0, 0)  # initialize documenting the current window to check

        for line in vcf:
            if line.startswith("##"):
                vcf_header = vcf_header + line
                continue
            elif line.startswith("#"):
                vcf_header = vcf_header + line
                anc_header = "\t".join(
                    [
                        line.strip("# ").split("\t", 9)[item]
                        for item in [0, 1, 2, 3, 4, 9]
                    ]
                )

                for i in range(num_ancs):  # Write header into each output vcf
                    files[f"vcf{i}"].write(vcf_header)
                    files[f"dos{i}"].write(anc_header)
                    files[f"ancdos{i}"].write(anc_header)

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

                for i in range(num_ancs):  # Write entries into each output file
                    output_lines[f"vcf{i}"] = vcf_out
                    output_lines[f"dos{i}"] = dos_anc_out
                    output_lines[f"ancdos{i}"] = dos_anc_out

                # optimized for quicker runtime - only move to next line when out of the current msp window
                # saves the current line until out of the window, then checks next line. VCF and window switches file should be in incremental order.
                while not (row[0] == window[0] and (window[1] <= pos < window[2])):
                    ancs = mspfile.readline()
                    if ancs.startswith("#"):  # skip the header lines
                        continue
                    if not ancs:
                        break  # when get to the end of the msp file, stop
                    # chm, spos, epos, sgpos, egpos, nsnps, calls
                    ancs_entry = ancs.strip().split("\t", 6)
                    calls = ancs_entry[6].split("\t")
                    window = (ancs_entry[0], int(ancs_entry[1]), int(ancs_entry[2]))

                # index by the number of individuals in the VCF file, should be half the number in the calls file
                for i, geno in enumerate(genos):
                    geno = geno.split(":")[0].split(
                        "|"
                    )  # assert incase eagle leaves some genos unphased
                    geno_a = str(geno[0])
                    geno_b = str(geno[1])
                    call_a = str(calls[2 * i])
                    call_b = str(calls[2 * i + 1])
                    counts = {anc: 0 for anc in range(num_ancs)}
                    anc_counts = {anc: 0 for anc in range(num_ancs)}
                    for j in range(num_ancs):
                        pop_genos[j] = ""
                        if call_a == str(j):
                            pop_genos[j] += geno_a
                            anc_counts[j] += 1
                            if geno_a == "1":
                                counts[j] += 1
                        else:
                            pop_genos[j] += "."

                        if call_b == str(j):
                            pop_genos[j] += "|" + geno_b
                            anc_counts[j] += 1
                            if geno_b == "1":
                                counts[j] += 1
                        else:
                            pop_genos[j] += "|."

                        output_lines[f"vcf{j}"] += "\t" + pop_genos[j]
                        output_lines[f"dos{j}"] += "\t" + str(counts[j])
                        output_lines[f"ancdos{j}"] += "\t" + str(anc_counts[j])

                for j in range(num_ancs):
                    output_lines[f"vcf{j}"] += "\n"
                    output_lines[f"dos{j}"] += "\n"
                    output_lines[f"ancdos{j}"] += "\n"

                    files[f"vcf{j}"].write(output_lines[f"vcf{j}"])
                    files[f"dos{j}"].write(output_lines[f"dos{j}"])
                    files[f"ancdos{j}"].write(output_lines[f"ancdos{j}"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--msp",
        help="path stem to RFmix msp file, not including .msp.tsv",
        required=True,
    )
    parser.add_argument(
        "--vcf-prefix",
        help="path stem to RFmix input VCF with phased genotypes, not including .vcf suffix",
        required=True,
    )
    parser.add_argument("--zipped", help="Input VCF is gzipped", action="store_true")
    parser.add_argument(
        "--num-ancs",
        type=int,
        help="Number of continental ancestries in admixed populations",
        default=2,
    )
    args = parser.parse_args()
    extract_tracts(**vars(args))
