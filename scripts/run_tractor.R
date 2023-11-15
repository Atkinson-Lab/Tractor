#!/usr/bin/env Rscript

# version: 1.1.0
script_version <- "1.1.0"
cat(paste("Tractor Script Version:", script_version), "\n")

# Changelog:
# Python script to R script to better deal with overly-inflated values.
# This script will instead output NA's in such cases.

require("optparse")

option_list = list(
  make_option(c("--hapdose"), type="character", default=NULL, 
              help="The prefix of hapcount and dosage files", metavar="character"),
  make_option(c("--phe"), type="character", default=NULL, 
              help="Phenotype and covariates file; 
                1st column sample ID, 2nd column phenotype, other columns will be treated as covariates; 
                The table should only contain numeric/integer values;
                Missing data is allowed", metavar="character"),
  make_option(c("--method"), type="character", default=NULL, 
              help="Either <linear> or <logistic>", metavar="character"),
  make_option(c("--out"), type="character", default=NULL, 
              help="Summary statistics file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



if (is.null(opt$hapdose)){
  print_help(opt_parser)
  stop("Missing hapdose files", call.=FALSE)
} else if (is.null(opt$phe)){
  print_help(opt_parser)
  stop("Missing phenotype/covariates file", call.=FALSE)
} else if (!(opt$method %in% c("linear", "logistic"))){
  print_help(opt_parser)
  stop("Method should be either <linear> or <logistic>", call.=FALSE)
} else if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Missing output file", call.=FALSE)
}


####### helper function
subset_mat_NA = function(rows, mat){
  t(sapply(rows, function(row, mat){if (row %in% row.names(mat)){return(mat[row,])}else{return(rep(NA, ncol(mat)))}},mat))
}

###########
RunTractor <- function(prefix, phefile, method, outfile){
  hapFiles = Sys.glob(paste0(prefix, ".*.hapcount.txt"))
  doseFiles = Sys.glob(paste0(prefix, ".*.dosage.txt"))
  inFiles = c(hapFiles, doseFiles)
  nAnc = length(doseFiles)
  phe = read.csv(phefile, sep = "\t", header = T)
  
  for (i in inFiles){
    cat("---------------------------\n")
    cat( paste("Load:", i, "\n"))
  }
  cat("Running Tractor... \n")
  
  ID = phe$IID 
  y = phe$y
  
  if (ncol(phe) > 2){
    COV_ = as.matrix(phe[,3:ncol(phe)])
  } else {
    COV_ = NULL
  }
  
  LAcolnames = paste0("LA", 0:(nAnc-1))
  Gcolnames = paste0("G", 0:(nAnc-1))

  
  data = lapply(inFiles,function(file){readLines(file)})
  nSNP = length(data[[1]]) - 1
  
  hap0ID = unlist(strsplit(data[[1]][1], "\t"))
  if (! all(ID == hap0ID[6:length(hap0ID)])){
    stop("sample in hapdose files doesn't match with the phenotype file", call.=FALSE)
  }
  
  
  resDF = setNames(data.frame(matrix(data = NA, nrow = 0, ncol = (5 +  4 * nAnc + 2 * (nAnc - 1)))),
                   c("CHR", "POS", "ID", "REF", "ALT",                                                                                             paste0("AF_anc",0:(nAnc-1)),
                     paste0("LAprop_anc",0:(nAnc-1)),
                     paste0("LAeff_anc",0:(nAnc-2)),
                     paste0("LApval_anc",0:(nAnc-2)),
                     paste0("Geff_anc",0:(nAnc-1)),
                     paste0("Gpval_anc",0:(nAnc-1))))
  write.table(resDF, outfile,  quote = F, row.names = F, sep = "\t")
  
  for (i in 1:nSNP){
    # matrix of Local Ancestry and Genotype
    LAG_ = matrix(unlist(lapply(data, function(file){
      fileSplit = unlist(strsplit(file[i+1], "\t"))
      as.integer(fileSplit[6:length(fileSplit)])
    })), nrow = nrow(phe), ncol = length(data))
    

    
    colnames(LAG_) = c(LAcolnames, Gcolnames) 
    
    AF = as.numeric(colSums(LAG_[,Gcolnames])/colSums(LAG_[,LAcolnames]))
    LAprop = as.numeric(colSums(LAG_[,LAcolnames])/sum(LAG_[,LAcolnames]))
    
    
    
    LAG_ = LAG_[,c(LAcolnames[1:(length(LAcolnames) - 1)], Gcolnames)]
    coef_rownames = paste0("LAG_", colnames(LAG_))
    LA_rownames = paste0("LAG_", LAcolnames[1:(length(LAcolnames) - 1)])
    G_rownames = paste0("LAG_", Gcolnames)
    
    if (method == "linear"){
        if (!is.null(COV_)){
          coefs = summary(lm(y ~ LAG_ + COV_))$coefficients
        } else {
          coefs = summary(lm(y ~ LAG_))$coefficients
        }
      
        reg_res = subset_mat_NA(coef_rownames, coefs[,c(1,4)])
        LAeff = as.numeric(reg_res[LA_rownames,1])
        LApval = as.numeric(reg_res[LA_rownames,2])
        Geff = as.numeric(reg_res[G_rownames,1])
        Gpval = as.numeric(reg_res[G_rownames,2])
        
        
    } else if (method == "logistic"){
        if (!is.null(COV_)){
          model = glm(y ~ LAG_ + COV_, family = binomial(link = "logit"))
          coefs = summary(model)$coefficients
        }else {
          model = glm(y ~ LAG_, family = binomial(link = "logit"))
          coefs = summary(model)$coefficients
        }
      
        reg_res = subset_mat_NA(coef_rownames, coefs[,c(1,4)])
        if (model$converged == TRUE){
          LAeff = as.numeric(reg_res[LA_rownames,1])
          LApval = as.numeric(reg_res[LA_rownames,2])
          Geff = as.numeric(reg_res[G_rownames,1])
          Gpval = as.numeric(reg_res[G_rownames,2])
        } else {
          # if glm doesn't converge, set effect size and P value as NA
          LAeff = rep(NA, length(LA_rownames))
          LApval = rep(NA, length(LA_rownames))
          Geff = rep(NA, length(G_rownames))
          Gpval = rep(NA, length(G_rownames))
        }
    }
    
    
    resDF[1, c("CHR", "POS", "ID", "REF", "ALT")] = unlist(strsplit(data[[1]][i+1], "\t"))[1:5]
    resDF[1, paste0("AF_anc",0:(nAnc-1))] = round(AF,6)
    resDF[1, paste0("LAprop_anc",0:(nAnc-1))] = round(LAprop,6)
    resDF[1, paste0("LAeff_anc",0:(nAnc-2))] = round(LAeff, 6)
    resDF[1, paste0("LApval_anc",0:(nAnc-2))] = LApval
    resDF[1, paste0("Geff_anc",0:(nAnc-1))] = round(Geff,6)
    resDF[1, paste0("Gpval_anc",0:(nAnc-1))] = Gpval
    
    write.table(resDF,  outfile, quote = F, row.names = F, col.names = F, append = T, sep = "\t")
  }
  
}

RunTractor(prefix = opt$hapdose, phefile = opt$phe, method = opt$method, outfile = opt$out)

