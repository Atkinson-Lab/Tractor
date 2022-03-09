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

def test_func(hapdose: str, phe: str, method: str, out:str):

    # read files 
    hapfiles = sorted(glob.glob(hapdose + '*.hapcount.txt'))
    dosefiles = sorted(glob.glob(hapdose + '*.dosage.txt'))

    file_list = hapfiles[:-1] + dosefiles
    print("Reading files....")
    for i in file_list:
        print("------")
        print(i)

    nParam = int((len(file_list) + 1)/2)

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

    



    if (method == "linear"):
        print("Runing Linear Regression Analysis...")
        with contextlib.ExitStack() as stack, open(out, 'w') as outfile:
            files = [stack.enter_context(open(fname)) for fname in file_list]
            
            # possible to add function to check errors
            lines = [file.readline().rstrip('\n').split(",") for file in files]
            
            if lines[0][0].split("\t")[5:] != list(phetbl.IID):
                print("ERROR: Phenotype ID must match with Hapdose file ID")
                return
            

            outfile.write(WriteHeader(nParam))  

            # use the first file as line to iterate
            while lines[0]:
                try:
                    lines = [file.readline() for file in files]
                    plines = [line.rstrip('\n').split("\t") for line in lines]
                    
                    variant = plines[0][0:5]
                    entry = np.array([Str2Int(line[5:]) for line in plines])

                    X1 = np.column_stack((np.transpose(entry), cov))

                    # use statsmodel at the moment, maybe will need to change it due to performance issue
                    X2 = sm.add_constant(X1)

                    # check singularity before running the program -- QR decomposition 
                    colboo = FindIndCol(X2)
                    X4model = X2[:,colboo]

                    # fit linear regression
                    est = sm.OLS(y, X4model, missing = 'drop').fit()


                    # two vectors that contains effect size estimates and p values
                    betaEst = np.full([len(colboo)], np.nan)
                    betaEst[colboo] = est.params

                    pEst = np.full([len(colboo)], np.nan)
                    pEst[colboo] = est.pvalues
                    betaEst[np.isnan(pEst)] = np.nan

                    # write betaEst and pEst into a new vector 
                    estimates = np.full([nParam*2], np.nan)
                    estimates[::2] = betaEst[nParam:2*nParam]
                    estimates[1::2] = pEst[nParam:2*nParam]
                    
                    outfile.write('\t'.join(plines[0][0:5]) +"\t" + "\t".join([str(int) for int in list(estimates)]) +"\n" )     
                    
                except:
                    print("END of calculation")


    elif (method == "logistic"):
        print("Runing Logistic Regression Analysis...")
        with contextlib.ExitStack() as stack, open(out, 'w') as outfile:
            files = [stack.enter_context(open(fname)) for fname in file_list]
            
            # possible to add function to check errors
            lines = [file.readline().rstrip('\n').split(",") for file in files]

            if lines[0][0].split("\t")[5:] != list(phetbl.IID):
                print("ERROR: Phenotype ID must match with Hapdose file ID")
                return

            outfile.write(WriteHeader(nParam))  

            # use the first file as line to iterate
            while lines[0]:
                try:
                    lines = [file.readline() for file in files]
                    plines = [line.rstrip('\n').split("\t") for line in lines]
                    

                    variant = plines[0][0:5]
                    entry = np.array([Str2Int(line[5:]) for line in plines])

                    X1 = np.column_stack((np.transpose(entry), cov))

                    try:
                        # use statsmodel at the moment, maybe will need to change it due to performance issue
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            X2 = sm.add_constant(X1)

                            # check singularity before running the program -- QR decomposition 
                            colboo = FindIndCol(X2)
                            X4model = X2[:,colboo]

                            est = sm.Logit(y, X4model, missing = 'drop').fit(disp=0)

                            # two vectors that contains effect size estimates and p values
                            betaEst = np.full([len(colboo)], np.nan)
                            betaEst[colboo] = est.params

                            pEst = np.full([len(colboo)], np.nan)
                            pEst[colboo] = est.pvalues
                            betaEst[np.isnan(pEst)] = np.nan

                        
                        estimates = np.full([nParam*2], np.nan)
                        estimates[::2] = betaEst[nParam:2*nParam]
                        estimates[1::2] = pEst[nParam:2*nParam]
                    
                        outfile.write('\t'.join(plines[0][0:5]) +"\t" + "\t".join([str(int) for int in list(estimates)]) +"\n" )  
                    except:
                        plines = [line.rstrip('\n').split("\t") for line in lines]
                        estimates = np.full([nParam*2], np.nan)
                        outfile.write('\t'.join(plines[0][0:5]) +"\t" + "\t".join([str(int) for int in list(estimates)]) +"\n" )  
                except:   
                    print("END of calculation")
                    

            


    else:
        print("Invalid --method input; Program terminating...")






# a function that convert str to int
def Str2Int(lst):
    return list(map(int, lst))

# a function to write header based on number of input files
def WriteHeader(nParam):
    header = "CHROM\tPOS\tID\tREF\tALT\t"
    for i in range(nParam):
        header = header + ("ANC{0}EFF\tANC{0}P\t".format(str(i))) 
    return(header[:-1] + "\n")

# a function that returns a boolean that represents which columns to include (check singularity)
def FindIndCol(X):
    return (np.invert(np.all(abs(np.linalg.qr(X)[1]) < 1e-10, axis=1)))


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
    test_func(**vars(args))
