% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_binary.R
\name{predict_binary}
\alias{predict_binary}
\title{predict_binary}
\usage{
predict_binary(
  motifs,
  xgb_model,
  celltype_prefix = NULL,
  pred = NULL,
  training_set = NULL,
  title = NULL
)
}
\arguments{
\item{motifs}{File containing motif counts}

\item{xgb_model}{File name of model}

\item{pred}{File name to output predicted values (default "predictions.txt")}

\item{training_set}{File name to output training counts (optional)}
}
\description{
Function to predict a set of test CREs using a binary classification model
}
\examples{
\dontrun{
extdata_path <- system.file("extdata",package = "BagOfMotifs")
motif_counts <- paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other_counts.txt")


predict_binary(motifs = motif_counts, xgb_model = paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other.rds"))
}
}
