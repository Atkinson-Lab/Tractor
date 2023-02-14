# Tractor script to run linear regression and logistic regression
__author__ = "TaotaoTan"


import argparse
import contextlib
import gzip
import logging
import glob
import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)



USAGE = """
RunTractor.py --hapdose  <prefix of deconvoluted ancestries and risk alleles files, XXXXX.ancN.hapcount.txt and XXXXX.ancN.dosage.txt>
                --phe <phenotype and covariates file; 1st column sample ID, 2nd column phenotype, other columns will be treated as covariates>
                --method <linear or logistic>
                --out <Output file>
"""

# user should make sure that the order of sample names matches

def run_tractor(hapdose: str, phe: str, method: str, out:str):

    # read files
    hapfiles = sorted(glob.glob(hapdose + '.*.hapcount.txt'))
    dosefiles = sorted(glob.glob(hapdose + '.*.dosage.txt'))
    file_list = hapfiles + dosefiles
    nParam = int(len(file_list)/2)
    print("Reading files....")
    for i in file_list:
        print("------")
        print(i)
        
    print("------")
    print("Notice:")
    print("Tractor drop one local ancestry term for regression. Therefore,", hapfiles[-1], "will not be used.")

    # read phenotype
    print("------")
    phetbl = pd.read_csv(phe, sep="\t")
    y = np.array(phetbl.iloc[:,1])
    cov = np.array(phetbl.iloc[:,2:])

    if (phetbl.columns[0] != "IID") | (phetbl.columns[1] != "y"):
        print("ERROR: Make sure include header your Phenotype file")
        print("Format should look like: ")
        print("IID\ty\t...")
        return
    
    if (method == "logistic") & (not(all(idx in (0.0, 1.0) for idx in phetbl.y[~np.isnan(phetbl.y)]))):
        print("ERROR: Phenotype must be encode as 0, 1, or blank")
        return

    with contextlib.ExitStack() as stack, open(out, 'w') as outfile:
        files = [stack.enter_context(open(fname)) for fname in file_list]

        # possible to add function to check errors
        lines = [file.readline().rstrip('\n').split(",") for file in files]

        if lines[0][0].split("\t")[5:] != list(phetbl.IID):
            print("ERROR: Phenotype ID must match with Hapdose file ID")

        outfile.write(WriteHeader(nParam))  

        # use the first file as line to iterate
        while lines[0]:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    lines = [file.readline() for file in files]
                    plines = [line.rstrip('\n').split("\t") for line in lines]
                    variant = plines[0][0:5]
                    entry = np.array([Str2Int(line[5:]) for line in plines])
                    AF_LAprop = Compute_AF_LAprop(entry, nParam)
                    X1 = np.column_stack((np.transpose(entry[np.r_[0:(nParam - 1), nParam:(2 * nParam)], ]), cov))
                    try:
                        X2 = sm.add_constant(X1)
                        
                        if (method == "linear"):
                            est = sm.OLS(y, X2, missing = 'drop').fit()
                        elif (method == "logistic"):
                            est = sm.GLM(y, X2, family=sm.families.Binomial(), missing = 'drop').fit(disp=0)
                        else:
                            print("ERROR: method must be either linear or logistic")
                            return


                        LAeff = est.params[1:nParam]
                        LApval = est.pvalues[1:nParam]
                        Geff = est.params[(nParam):(2*nParam)]
                        Gpval = est.pvalues[(nParam):(2*nParam)]
                        
                        # fix coefficients if allele frequency = 0
                        coef2fix_af0 = abs(np.array(AF_LAprop[0:nParam])) < 1e-10
                        if np.any(coef2fix_af0):
                            Geff[coef2fix_af0] = np.nan
                            Gpval[coef2fix_af0] = np.nan
                        
                        # fix coefficients if allele frequency = 1
                        # sometimes the program will explode
                        coef2fix_af1 = abs(np.array(AF_LAprop[0:nParam]) - 1 ) < 1e-10
                        if (coef2fix_af1[nParam - 1] == True) or (np.sum(coef2fix_af1) > 1) or np.isnan(est.llf):
                            Geff[0:nParam] = np.nan
                            Gpval[0:nParam] = np.nan
                            LAeff[0:(nParam-1)] = np.nan
                            LApval[0:(nParam-1)] = np.nan
                        else:
                            Geff[coef2fix_af1] = Geff[coef2fix_af1]*2
                            LAeff[coef2fix_af1[0:(nParam-1)]] = np.nan
                            LApval[coef2fix_af1[0:(nParam-1)]] = np.nan
                        outfile.write('\t'.join(plines[0][0:5]) +"\t" + 
                                    "\t".join([str(k) for k in AF_LAprop]) + "\t" +
                                    "\t".join([str(k) for k in list(LAeff)]) +"\t" + 
                                    "\t".join([str(k) for k in list(LApval)]) + "\t" +
                                    "\t".join([str(k) for k in list(Geff)]) +"\t" + 
                                    "\t".join([str(k) for k in list(Gpval)]) + "\n" )  
                        
                    except:
                        lines = [file.readline() for file in files]
                        plines = [line.rstrip('\n').split("\t") for line in lines]
                        variant = plines[0][0:5]
                        entry = np.array([Str2Int(line[5:]) for line in plines])
                        AF_LAprop = Compute_AF_LAprop(entry, nParam)
                        Geff = np.full([nParam], np.nan)
                        Gpval = np.full([nParam], np.nan)
                        LAeff = np.full([nParam - 1], np.nan)
                        LApval = np.full([nParam - 1], np.nan)
                        outfile.write('\t'.join(plines[0][0:5]) +"\t" + 
                        "\t".join([str(k) for k in AF_LAprop]) + "\t" +
                        "\t".join([str(k) for k in list(LAeff)]) +"\t" + 
                        "\t".join([str(k) for k in list(LApval)]) + "\t" +
                        "\t".join([str(k) for k in list(Geff)]) +"\t" + 
                        "\t".join([str(k) for k in list(Gpval)]) + "\n" )  

            except:   
                print("END of calculation")





# a function to write the header
def WriteHeader(nParam):
    terms = ["AF", "LAprop", "LAeff", "LApval", "Geff", "Gpval"]
    append1 = ["_anc"+ str(k) for k in list(range(0, nParam))]
    append2 = ["_anc"+ str(k) for k in list(range(0, nParam - 1))]
    header = "CHROM\tPOS\tID\tREF\tALT\t"
    for i in terms:
        if (i == "LAeff" or i == "LApval"):
            header = header + '\t'.join([i + j for j in append2]) + "\t"
        else:
            header = header + '\t'.join([i + j for j in append1]) + "\t"
    header = header + "\n"
    return header

def Str2Int(lst):
    return list(map(int, lst))

def Compute_AF_LAprop(entry, nParam):
    LAcount = np.sum(entry[0:nParam, :], axis = 1)
    Gtot = np.sum(entry[nParam:(2*nParam + 1), :], axis = 1)
    LAtot = sum(LAcount)
    LAprop = LAcount/LAtot
    AF = Gtot/LAcount
    return (AF.round(5).tolist() + LAprop.round(5).tolist() )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hapdose",
        help="path stem to deconvoluted ancestries and risk alleles files, not including .ancN.hapcount.txt or .ancN.dosage.txt",
        required=True,
    )
    parser.add_argument(
        "--phe",
        help="phenotype and covariate file; 1st column sampleID, 2nd column phenotype, rest covariates",
        required=True,
    )
    parser.add_argument("--method", help="Whether perform linear regression or logistic regression", required=True)
    parser.add_argument("--out", help="output file for Tractor sumstats", required=True)

    args = parser.parse_args()
    run_tractor(**vars(args))
