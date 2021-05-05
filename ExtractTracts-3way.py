# Tractor script to extract out ancestry segments from each reference population from admixed samples.
__author__ = "egatkinson"


import argparse
import gzip
import contextlib

USAGE = """
ExtractTracts.py --msp  <an ancestral calls file produced by RFmix version 2, suffixed with .msp.tsv>
                             --vcf <VCF file suffixed with .vcf>
"""

# input is expected to be a VCF file suffixed with .vcf and an ancestral calls file produced by RFmix version 2, suffixed with .msp.tsv
# output will be returned in the same location as input VCF file with same naming convention


def extract_tracts(msp, vcf_prefix, zipped, num_ancs):
    mspfile = f"{msp}.msp.tsv"
    vcf = f"{vcf_prefix}.vcf.gz" if zipped else f"{vcf_prefix}.vcf"
    output_files = {}

    for i in range(num_ancs):
        output_files[
            f"out{i}"
        ] = f"{vcf_prefix}.anc{i}.vcf"  # output extracted VCF for each ancestry into separate files
        output_files[
            f"outdos{i}"
        ] = f"{vcf_prefix}.anc{i}.dosage.txt"  # output dosages for each ancestry into separate files
        output_files[
            f"outancdos{i}"
        ] = f"{vcf_prefix}.anc{i}.hapcount.txt"  # output number of haplotype for each ancestry into separate files

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
                    files[f"out{i}"].write(vcf_header)
                    files[f"outdos{i}"].write(anc_header)
                    files[f"outancdos{i}"].write(anc_header)

            if not line.startswith("#"):
                # Entry format is ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'genotypes']
                row = line.strip().split("\t", 9)

                # print(row)
                # Grab fields needed for each file output
                vcf_out = "\t".join(row[:9])
                dos_anc_out = "\t".join(row[:5])

                # split each strand geno call apart from each other
                genos = row[9].split("\t")
                pos = int(row[1])
                output_lines = {}
                pop_genos = {}

                for i in range(num_ancs):  # Write entries into each output file
                    output_lines[f"output{i}"] = vcf_out
                    output_lines[f"outputdos{i}"] = dos_anc_out
                    output_lines[f"outputancdos{i}"] = dos_anc_out
                    pop_genos[i] = ""

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

                for i, geno in enumerate(
                    genos
                ):  # index by the number of individuals in the VCF file, should be the same number in the calls file
                    geno = geno.split(":")[0].split("|")
                    geno_a = str(geno[0])
                    geno_b = str(geno[1])
                    call_a = str(calls[2 * i])
                    call_b = str(calls[2 * i + 1])
                    counts = {anc: 0 for anc in range(num_ancs)}
                    anc_counts = {anc: 0 for anc in range(num_ancs)}
                    for j in range(num_ancs):
                        pop_genos[j] = ""
                        if call_a == str(j):  # j may need to be converted to string
                            pop_genos[j] += geno_a
                            anc_counts[j] += 1
                            if geno_a == "1":
                                counts[j] += 1
                        else:
                            pop_genos[j] += "."

                        if call_b == str(j):  # j may need to be converted to string
                            pop_genos[j] += "|" + geno_b
                            anc_counts[j] += 1
                            if geno_b == "1":
                                counts[j] += 1
                        else:
                            pop_genos[j] += "|."

                        output_lines[f"output{j}"] += "\t" + pop_genos[j]
                        output_lines[f"outputdos{j}"] += "\t" + str(counts[j])
                        output_lines[f"outputancdos{j}"] += "\t" + str(anc_counts[j])

                for j in range(num_ancs):
                    output_lines[f"output{j}"] += "\n"
                    output_lines[f"outputdos{j}"] += "\n"
                    output_lines[f"outputancdos{j}"] += "\n"

                    files[f"out{j}"].write(output_lines[f"output{j}"])
                    files[f"outdos{j}"].write(output_lines[f"outputdos{j}"])
                    files[f"outancdos{j}"].write(output_lines[f"outputancdos{j}"])


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
        help="Number of continental ancestries in admixed populations",
        default=2,
    )
    args = parser.parse_args()
    extract_tracts(**vars(args))
