#!/usr/bin/env Rscript

# version: 1.4.0
script_version <- "1.4.0"
cat(paste("Tractor Script Version :", script_version), "\n")

# Changelog

## v1.4.0  -  2024-05-10 (Enhanced File Handling, Performance Improvement)
### - Added support for compressed (gz) hapcount/dosage and phenotype files.
### - Improved file reading efficiency by implementing fread in chunks, mitigating memory errors.
### - Implemented parallel processing for regression using doParallel, resulting in significant speed improvements with multi-core systems.
### - Enhanced flexibility in organizing phenotype files:
###     - Users can specify sample ID column, phenotype ID column, and covariate column list.
### - Updated output summary statistics to include SE and t-val, with column names adjusted to adhere to GWAS standards.

## v1.2.0 and v1.3.0 - Internal Updates

## v1.1.0 2024-02-09 (Bug Fixes and Script Conversion)
### - Fixed an issue with inflated values from v1.0.0 by converting the Python script to R
### - Modified the script to output NA's instead of erroneous values in such cases.

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressWarnings(suppressMessages(library(R.utils)))
suppressMessages(library(dplyr))
suppressMessages(library(doParallel))

option_list = list(
  make_option(c("--hapdose"), type="character", default=NULL, metavar="character",
              help="[Mandatory] Prefix for hapcount and dosage files.
                E.g. If you have the following files: filename.anc0.dosage.txt, filename.anc0.hapcount.txt, filename.anc1.dosage.txt, filename.anc1.hapcount.txt,
                use \"--hapdose filename\""),
  make_option(c("--phenofile"), type="character", default=NULL, metavar="character",
              help="[Mandatory] Specify the file path containing phenotype and covariate data.
                Default assumptions: Sample ID column: \"IID\" or \"#IID\", Phenotype column: \"y\".
                If different column names are used, refer to --sampleidcol and --phenocol arguments.
                All covariates MUST be included using --covarcollist."),
  make_option(c("--covarcollist"), type="character",metavar="character",
              help="[Mandatory] Specify column names of covariates in the --phenofile.
                Only listed columns will be included as covariates. Separate multiple covariates with commas.
                E.g. --covarcollist age,sex,PC1,PC2.
                To exclude covariates, specify \"--covarcollist none\"."),
  make_option(c("--method"), type="character", default=NULL, metavar="character",
              help="[Mandatory] Specify the method to be used: <linear> or <logistic>."),
  make_option(c("--output"), type="character", default=NULL, metavar="character",
              help="[Mandatory] File name for summary statistics output.
                E.g. /path/to/file/output_sumstats.txt"),
  make_option(c("--sampleidcol"), type="character", default=NULL,metavar="character",
              help="[Optional] Specify sample ID column name in the --phenofile.
                Default: \"IID\" or \"#IID\""),
  make_option(c("--phenocol"), type="character", default=NULL,metavar="character",
              help="[Optional] Specify phenotype column name in the --phenofile.
                Default: \"y\""),
  make_option(c("--chunksize"), type="numeric", default=10000, metavar="numeric",
              help="[Optional] Number of rows to read at once from hapcount and dosage files.
                Use smaller values for lower memory usage.
                Default: 10000
                Note: Higher chunksize speeds up streaming but requires more memory.
                      If out-of-memory errors occur, try reducing --chunksize or --nthreads (if unable to increase memory)."),
  make_option(c("--nthreads"), type="numeric", default=1, metavar="numeric",
              help="[Optional] Specify number of threads to use.
                Increasing threads can speed up processing but may increase memory usage.
                Default: 1"),
  make_option(c("--totallines"), type="numeric", default=NULL, metavar="numeric",
              help="[Optional] Specify total number of lines in hapcount/dosage files (can be verified using wc -l).
                If not provided, it will be calculated internally (recommended).
                Exercise caution: if --totallines is smaller than the actual lines in the hapcount/dosage files, 
                only a subset of data will be analyzed. If larger than the actual lines in the files,
                an error will occur. Both scenarios are discouraged.")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$hapdose)) {
  print_help(opt_parser)
  stop("Error: Missing prefix for --hapdose files. Check --help", call.=FALSE)
} else if (!(opt$method %in% c("linear", "logistic"))) {
  print_help(opt_parser)
  stop("Error: Method should be either <linear> or <logistic>", call.=FALSE)
} else if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Error: Missing output file (the file in which summary statistics will be saved). E.g. <prefix>_sumstats.txt", call.=FALSE)
} else if (is.null(opt$phenofile) || !file.exists(opt$phenofile)) {
  print_help(opt_parser)
  stop("Error: Missing or can't locate phenotype file. Check --phenofile argument again.")
} else if (is.null(opt$covarcollist)) {
  print_help(opt_parser)
  stop("Error: Missing values for covariates. Check --covarcollist (mandatory) argument. If you wish to use no covariates, specify \"--covarcollist none\".")
}

# Helper functions
# Get necessary values from regression summary
subset_mat_NA = function(rows, mat) {
  rows %>%
    sapply(function(row, mat) {
      if (row %in% row.names(mat)) {
        return(mat[row, ])
      } else {
        return(rep(NA, ncol(mat)))
      }
    }, mat = mat) %>%
    t()
}

# Extract model coefficients and values
extract_model_info <- function(model, coef_rownames, LA_rownames, G_rownames) {
  
  coefs     <- summary(model)$coefficients
  reg_res   <- subset_mat_NA(coef_rownames, coefs)
  
  LAeff     <- as.numeric(reg_res[LA_rownames, 1])
  LApval    <- as.numeric(reg_res[LA_rownames, 4])
  
  Geff      <- as.numeric(reg_res[G_rownames, 1])
  Gstderr   <- as.numeric(reg_res[G_rownames, 2])
  Gtval     <- as.numeric(reg_res[G_rownames, 3])
  Gpval     <- as.numeric(reg_res[G_rownames, 4])
  
  N         <- sum(complete.cases(model$model))
  
  # Return a list with extracted information
  return(list(
    LAeff  = LAeff,
    LApval = LApval,
    Geff   = Geff,
    Gstderr = Gstderr,
    Gtval  = Gtval,
    Gpval  = Gpval,
    N      = N
  ))
}

RunTractor <- function(prefix, phenofile, sampleidcol, phenocol, covarcollist, chunksize, totallines, method, outfile, nthreads) {
  
  hapFiles = Sys.glob(paste0(prefix, ".*.hapcount.txt*"))
  doseFiles = Sys.glob(paste0(prefix, ".*.dosage.txt*"))
  if (length(hapFiles)!=length(doseFiles)) {
    cat("Hapcount Files         : ", hapFiles, "\n")
    cat("Dosage Files           : ", doseFiles, "\n")
    stop("Error: Reading different number of hapcount and dosage files. The number of hapcount and dosage files should be identical.")
  }
  inFiles = c(hapFiles, doseFiles)
  txt_files <- sum(endsWith(inFiles, ".txt"))
  gz_files  <- sum(endsWith(inFiles, ".gz"))
  
  if (txt_files > 0 && gz_files > 0) {
    stop("Error: Both .txt and .gz hapcount/dosage files were detected, which is not recommended. Please make sure all the files are of the same format.")
  } else if (txt_files == 0 && gz_files == 0) {
    stop("Error: No hapcount/dosage files (txt or txt.gz) files were detected. Please make sure all necessary files are present in the directory.")
  }
  rm(txt_files)
  rm(gz_files)
  nAnc = length(doseFiles)
  df_phe = fread(phenofile, sep="\t", header=T)
  
  # sampleidcol
  if (!is.null(sampleidcol)) {
    if (sampleidcol %in% colnames(df_phe)) {
      colnames(df_phe)[colnames(df_phe) == sampleidcol] <- "IID"
      cat("Sample ID column used  : ", sampleidcol, "\n")
    } else {
      stop(paste("Error: Column", sampleidcol, "does not exist in the file."))
    }
  } else {
    if ("#IID" %in% colnames(df_phe)) {
      colnames(df_phe)[colnames(df_phe) == "#IID"] <- "IID"
      cat("Sample ID column used  : #IID\n")
    } else if (!("IID" %in% colnames(df_phe))) {
      stop(paste("Error: Unable to identify sample ID column. Default column name expected for sample ID is IID or #IID.
                  Alternatively, sample ID can be provided with --sampleidcol"))
    } else {
      cat("Sample ID column used  : IID\n")
    }
  }
  
  # phenocol
  if (!is.null(phenocol)) {
    if (phenocol %in% names(df_phe)) {
      colnames(df_phe)[colnames(df_phe) == phenocol] <- "y"
      cat("Phenotype column used  : ", phenocol, "\n")
    } else {
      stop("Error: ", phenocol, " not found in --phenofile. Please check column names again.")
    }
  } else if ("y" %in% names(df_phe)) {
    cat("Phenotype column used  : y\n")
  } else {
    stop(paste("Error: Unable to identify phenotype column. Default column name expected for phenotype is y.
                Alternatively, phenotype column name can be provided with --phenocol"))
  }
  
  # covarcollist
  if ( tolower(covarcollist) != "none") {
    covars <- strsplit(covarcollist,",")[[1]]
    for (i in c(covars)) {
      if (!(i %in% names(df_phe))) {
        stop("Error: Covariate column: ", i, " does not exist in the --phenofile.\n")
      }
    }
    cat("Covariates used        :", covars, "\n")
  } else {
    cat("Covariates used        : NO COVARIATES USED. If you wish to use covariates, please use the --covarcollist argument.\n")
  }
  
  cat("Total hapcount/dosage files identified  :", length(inFiles), "\n")
  for (i in inFiles){
    cat("                        ",i,"\n")
  }
  
  if (is.null(totallines)) {
    if (!endsWith(inFiles[1],".gz")) {
      wc_op <- paste("wc -l", inFiles[1]) %>%
        system(intern = T) %>%
        trimws("left") %>% 
        strsplit(" ")
    } else {
      wc_op <- paste("zcat <", inFiles[1], "| wc -l") %>%
        system(intern = T) %>%
        trimws("left") %>% 
        strsplit(" ")
    }
    
    totallines <- wc_op[[1]][1] %>% as.numeric() -1
    cat("Total SNPs Identified  :", totallines, "\n")
    rm(wc_op)
  }
  
  LAcolnames = paste0("LA", 0:(nAnc-1))
  Gcolnames = paste0("G", 0:(nAnc-1))
  COV_      = NULL
  
  iters          <- 1
  skip_val       <- 1
  max_iters      <- ceiling(totallines/chunksize)
  data_colnames  <- NULL
  
  # Set parallelization (Use nthreads-1)
  if(nthreads>detectCores()) {
    stop(paste0("Error: --nthreads must be less than ", detectCores(),". Provided --nthreads ", nthreads))
  }
  if (nthreads > 1) {
    nthreads = nthreads - 1
  } else if (nthreads <= 0) {
    stop(paste0("Error: --nthreads must be 1 or more. Provided --nthreads ", nthreads))
  }
  cl<-makeCluster(nthreads)
  registerDoParallel(cl)
  cat("Threads being used     :", getDoParWorkers(),"\n") # Print number of cores being used.
  
  dopar_packages <- c("data.table","dplyr")
  dopar_functions <- c("subset_mat_NA","extract_model_info")
  
  while (iters <= max_iters) {
    
    if (iters != 1) {
      data = lapply(inFiles,function(file) {
        data.table::fread(file, nrows=chunksize, skip=skip_val,
                          col.names=data_colnames, sep="\t") #, header=TRUE
      })
    } else {
      data = lapply(inFiles,function(file) {
        data.table::fread(file, nrows=chunksize, skip=skip_val-1,
                          sep="\t", header=TRUE)
      })
      data_colnames <- colnames(data[[1]])
      
      # Match phenotype sample ID order to hapdos sample ID order.
      sampleID_hapdos <- data.table(IID=colnames(data[[1]])[-(1:5)])
      
      df_phe$IID <- as.character(df_phe$IID)
      
      # Samples present in phenotype file but not in hapcount/dosage files
      excluded_samples <- anti_join(df_phe, sampleID_hapdos) %>% select("IID")
      path_to_excluded_samples <- file.path(dirname(outfile),"samples_excluded_from_phenotype.txt")
      if (nrow(excluded_samples) > 0) {
        fwrite(excluded_samples, path_to_excluded_samples,
               sep = "\t", quote = FALSE, col.names = FALSE)
        cat("NOTE: Some samples found in phenotype file were not present in hapcount and dosage files.\n")
        cat("      These samples are being excluded. For more information, check", path_to_excluded_samples,"\n")
      }
      
      # Left join to rearrange phenotype file to match order of sample IDs in hapcount/dosage files
      df_phe <- left_join(sampleID_hapdos, df_phe)
      ID <- df_phe$IID
      y <- df_phe$y
      
      cat("Total samples          :", nrow(df_phe),"(from hapcount/dosage files)\n")
      # Confirm order of ID matches that of data[[1]]
      if (! all(ID == colnames(data[[1]])[6:length(data[[1]])])){
        stop("Error: Samples in hapdose files doesn't match with the phenotype file", call.=FALSE)
      }
      rm(ID)
      
      # Create covariate matrix from phenotype file.
      if (tolower(covarcollist) != "none") {
        COV_ <- as.matrix(df_phe[, ..covars])
      } else {
        COV_ = NULL
      }
      rm(sampleID_hapdos)
      rm(excluded_samples)
      rm(covarcollist)
    }
    
    nSNP = nrow(data[[1]])
    results <- foreach (i=1:nSNP,
                        .combine='c',
                        .inorder = TRUE,
                        .packages=dopar_packages,
                        .export=dopar_functions) %dopar% {
                          
                          LAG_ = sapply(data, function(dt) dt[i, 6:ncol(dt), with = FALSE]) %>%
                            as.numeric() %>%
                            matrix(nrow=nrow(df_phe), ncol=length(data))
                          
                          colnames(LAG_) = c(LAcolnames, Gcolnames)
                          
                          AF            = as.numeric(colSums(LAG_[,Gcolnames])/colSums(LAG_[,LAcolnames]))
                          LAprop        = as.numeric(colSums(LAG_[,LAcolnames])/sum(LAG_[,LAcolnames]))
                          
                          # Sum of Local Ancestry equals 1, so exclude one LA term
                          LAG_          = LAG_[,c(LAcolnames[1:(length(LAcolnames) - 1)], Gcolnames)]
                          coef_rownames = paste0("LAG_", colnames(LAG_))
                          LA_rownames   = paste0("LAG_", LAcolnames[1:(length(LAcolnames) - 1)])
                          G_rownames    = paste0("LAG_", Gcolnames)
                          
                          if (method == "linear") {
                            if (!is.null(COV_)){
                              model = lm(y ~ LAG_ + COV_)
                            } else {
                              model = lm(y ~ LAG_)
                            }
                            
                            extracted_info <- extract_model_info(model, coef_rownames, LA_rownames, G_rownames)
                            LAeff  = extracted_info$LAeff
                            LApval = extracted_info$LApval
                            Geff   = extracted_info$Geff
                            Gstderr = extracted_info$Gstderr
                            Gtval  = extracted_info$Gtval
                            Gpval  = extracted_info$Gpval
                            N      = extracted_info$N
                          } else if (method == "logistic") {
                            if (!is.null(COV_)){
                              model = glm(y ~ LAG_ + COV_, family = binomial(link = "logit"))
                            } else {
                              model = glm(y ~ LAG_, family = binomial(link = "logit"))
                            }
                            
                            if (model$converged == TRUE){
                              extracted_info <- extract_model_info(model, coef_rownames, LA_rownames, G_rownames)
                              LAeff  = extracted_info$LAeff
                              LApval = extracted_info$LApval
                              Geff   = extracted_info$Geff
                              Gstderr = extracted_info$Gstderr
                              Gtval  = extracted_info$Gtval
                              Gpval  = extracted_info$Gpval
                              N      = extracted_info$N
                            } else {
                              # if glm doesn't converge, set effect size and P value as NA
                              LAeff     = rep(NA, length(LA_rownames))
                              LApval    = rep(NA, length(LA_rownames))
                              Geff      = rep(NA, length(G_rownames))
                              Gstderr   = rep(NA, length(G_rownames))
                              Gtval     = rep(NA, length(G_rownames))
                              Gpval     = rep(NA, length(G_rownames))
                              N         = sum(complete.cases(model$model))
                            }
                          }
                          
                          temp_result <- data.table(CHR = data[[1]][[i,1]],
                                                    POS = data[[1]][[i,2]],
                                                    ID  = data[[1]][[i,3]],
                                                    REF = data[[1]][[i,4]],
                                                    ALT = data[[1]][[i,5]],
                                                    N   = N)
                          for (ancid in 0:(nAnc-1)) {
                            temp_result[[paste0("AF_anc",ancid)]]         <- round(AF[ancid+1],6)
                            temp_result[[paste0("LAprop_anc", ancid)]]    <- round(LAprop[ancid+1],6)
                            temp_result[[paste0("beta_anc", ancid)]]      <- round(Geff[ancid+1],6)
                            temp_result[[paste0("se_anc", ancid)]]        <- round(Gstderr[ancid+1], 4)
                            temp_result[[paste0("pval_anc", ancid)]]      <- Gpval[ancid+1]
                            temp_result[[paste0("tval_anc", ancid)]]      <- round(Gtval[ancid+1], 4)
                          }
                          for (ancid in 0:(nAnc-2)) {
                            temp_result[[paste0("LApval_anc",ancid)]]     <- LApval[ancid+1]
                            temp_result[[paste0("LAeff_anc", ancid)]]     <- round(LAeff[ancid+1], 6)
                          }
                          
                          # df_result <- rbind(df_result, temp_result)
                          list(temp_result)
                        }
    
    combined_result <- rbindlist(results)
    
    if (iters == max_iters) {
      # Last iteration - Append and stop threading/parallelization
      fwrite(combined_result, file = outfile, quote = F, append = T, sep = "\t")
      rm(combined_result)
      cat(paste0("Chunk ",iters,"/",max_iters, " completed.     ", Sys.time()),"\n")
      cat("Run complete.\n\n\n")
      stopCluster(cl) # Stop threading/parallelization
      break # End loop
    } else if (iters != 1) {
      # Append results to a pre-existing file
      fwrite(combined_result, file = outfile, quote = F, append = T, sep = "\t")
      rm(combined_result)
    } else {
      # Create, or overwrite pre-existing file in first iteration
      fwrite(combined_result, file = outfile, quote = F, sep = "\t")
      rm(combined_result)
      cat("\n")
    }
    
    cat(paste0("Chunk ",iters,"/",max_iters, " completed.     ", Sys.time()),"\n")
    
    # Updating looping variables
    iters    = iters + 1
    skip_val = skip_val + chunksize
  }
}

RunTractor(prefix = opt$hapdose,
           phenofile = opt$phenofile,
           sampleidcol = opt$sampleidcol,
           phenocol = opt$phenocol,
           covarcollist = opt$covarcollist,
           chunksize = opt$chunksize,
           totallines = opt$totallines,
           method = opt$method,
           outfile = opt$output,
           nthreads = opt$nthreads)

