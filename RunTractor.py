# Tractor script to run linear regression and logistic regression
__author__ = "TaotaoTan"


import argparse
import contextlib
import gzip
import logging
import glob
import numpy as np
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import pandas as pd
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

    if (method == "linear"):
        print("Runing Linear Regression Analysis...")
        with contextlib.ExitStack() as stack, open(out, 'w') as outfile:
            files = [stack.enter_context(open(fname)) for fname in file_list]
            
            # possible to add function to check errors
            lines = [file.readline().rstrip('\n').split(",") for file in files]

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
                    est = sm.OLS(y, X2).fit()
                    
                    estimates = np.full([nParam*2], np.nan)
                    estimates[::2] = est.params[nParam:2*nParam]
                    estimates[1::2] = est.pvalues[nParam:2*nParam]
                    
                    outfile.write('\t'.join(plines[0][0:5]) +"\t" + "\t".join([str(int) for int in list(estimates)]) +"\n" )     
                    
                except:
                    print("END of calculation")


    elif (method == "logistic"):
        print("Runing Logistic Regression Analysis...")
        with contextlib.ExitStack() as stack, open(out, 'w') as outfile:
            files = [stack.enter_context(open(fname)) for fname in file_list]
            
            # possible to add function to check errors
            lines = [file.readline().rstrip('\n').split(",") for file in files]

            outfile.write(WriteHeader(nParam))  

            # use the first file as line to iterate
            i = 1
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
                            est = sm.Logit(y, X2).fit(disp=0)
                        
                        estimates = np.full([nParam*2], np.nan)
                        estimates[::2] = est.params[nParam:2*nParam]
                        estimates[1::2] = est.pvalues[nParam:2*nParam]
                    
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
        header = header + ("ANC{0}EFF\tANCP{0}\t".format(str(i))) 
    return(header[:-1] + "\n")



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