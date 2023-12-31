% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_for_binary_model.R
\name{binModel}
\alias{binModel}
\title{binModel}
\usage{
binModel(
  data_path,
  qval_thresh,
  outDir,
  target_ct = NULL,
  nthreads = 1,
  nrounds = 10000,
  eta = 0.01,
  max_depth = 6,
  subsample = 0.5,
  colsample_bytree = 0.5,
  objective = "binary:logistic",
  early_stopping_rounds = NULL,
  eval_metric = "error",
  maximize = F,
  params = list(),
  feval = NULL,
  verbose = 0,
  print_every_n = 1L,
  save_period = NULL,
  xgb_model = NULL,
  callbacks = list(),
  training = 0.6
)
}
\arguments{
\item{data_path}{Path to the input data directory.}

\item{qval_thresh}{q-value threshold for motif filtering. Default to 0.5 (q-value <= 0.5).}

\item{outDir}{Name of directory to save output files. Output files will be named cellType_vsOthers.}

\item{target_ct}{Name of the target cell type/condition. If target_ct is NULL, the function will produce a matrix for each cell type in the input directory.}

\item{nthreads}{The number of threads to use}
}
\description{
Function to produce a matrix of motif frequency. Motif instances will be filtered using the q-value threshold provided. 
The regions annotated to the target cell type/condition will be coded as 1. Background regions will be coded as 0.
The outpud matrix will contain a balanced number of positive (target) and negative (background) instances.
}
\examples{
\dontrun{
extdata_path <- system.file("extdata",package = "BagOfMotifs")
data_path <- paste0(extdata_path, "/tutorial/motifs") 

binModel(data_path=data_path, qval_thresh=0.5, outDir='results/', target_ct= 'Cardiomyocytes')
}
}
