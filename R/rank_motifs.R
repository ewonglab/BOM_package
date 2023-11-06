#############################################################
#' rankMotifs
#' @description Function to rank motifs by the sum or mean of absolute SHAP scores.  
#' 
#' @param shap_file  Path to the file containing SHAP scores
#' @param out_file File name to output the motifs ranked by SHAP (default to "ranked_motifs")
#' @param rank_type Rank type. Either 'sum' or 'mean'
#' 
#' @examples 
#'
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' binPredictions <- paste0(extdata_path, "/tutorial/Cardiomyocytes_vs_other_SHAP.txt")
#' 
#' 
#' rankMotifs(shap_file = binPredictions, out_file = "ranked_motifs", type = "sum")
#' }
#' @export
rankMotifs <- function(shap_file, type, out_file){
  
    if(!type %in% c("sum", "mean")){
    stop(paste("Invalid ranking type:", type))
    }
  message("Reading SHAP values...\n")
  shap <- read.table(file = shap_file, header = T, stringsAsFactors = F, sep ='\t', row.names = 1)
  
  if(type == "sum"){
    shap_per_motif <- apply(shap, 2, function(x) sum(abs(x)))
  }
  if(type == "mean"){
    shap_per_motif <- apply(shap, 2, function(x) mean(abs(x)))
  }
  
  shap_per_motif <- shap_per_motif[order(-shap_per_motif)]
  shap_per_motif <- as.data.frame(shap_per_motif)
  
  if(type == "sum"){
    colnames(shap_per_motif) <- "sum_abs_SHAP"
  }
  if(type == "mean"){
    colnames(shap_per_motif) <- "mean_abs_SHAP"
  }
  
  # return(shap_per_motif)
  message("Saving ranked motifs...\n")
  write.table(x = shap_per_motif, file = out_file, quote = F)
  
}
