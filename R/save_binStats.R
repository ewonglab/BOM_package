#############################################################
## binStats_single
#' @description This function calculates prediction statistics for a single predictions file
#'
#' @param pred_df input data frame containing predictions from a binary model.
#' @param celltype Cell type/state.
#'
#' @examples
#'
#' Endo_pred <- read.table(file = "Endothelium_pred.txt", header = TRUE, stringsAsFactors = FALSE)
#' Endo_stats <- binStats_single(pred_df = Endo_pred, celltype = "Endothelium")
#'
#' @export
#'

binStats_single <- function(pred_df, model){
  TP <- nrow(pred_df[pred_df$true_class == 1 & pred_df$predicted_class == 1, ])
  TN <- nrow(pred_df[pred_df$true_class == 0 & pred_df$predicted_class == 0, ])
  FP <- nrow(pred_df[pred_df$true_class == 0 & pred_df$predicted_class == 1, ])
  FN <- nrow(pred_df[pred_df$true_class == 1 & pred_df$predicted_class == 0, ])
  recall <- TP/(TP + FN)
  precision <- TP/(TP + FP)
  f1 <- 2 * ((precision * recall)/(precision + recall))
  
  roc_auc <- round(cvAUC::AUC(pred_df$prob, pred_df$true_class), 4)
  acc <- round(nrow(pred_df[pred_df$predicted_class == pred_df$true_class, ]) / nrow(pred_df), 4)
  
  mcc <- ((as.numeric(TP) * as.numeric(TN)) - (as.numeric(FP) * as.numeric(FN))) /
    sqrt((as.numeric(TP) + as.numeric(FP)) * (as.numeric(TP) + as.numeric(FN)) * 
         (as.numeric(TN) + as.numeric(FP)) * (as.numeric(TN) + as.numeric(FN)))
  
  pred_df$true_class <- factor(pred_df$true_class, levels = c(0, 1))
  pr_auc_val <- yardstick::pr_auc(pred_df, truth = true_class
                                  , prob, event_level = "second"
                                  , estimator = "binary")
  pr_auc_val <- pr_auc_val$.estimate
  
  # output dataframe
  out_df <- data.frame(Model = model, Accuracy = acc
                       , auPR = pr_auc_val, auROC = roc_auc
                       , F1 = f1, MCC = mcc, Precision = precision
                       , Recall = recall)
  return(out_df)
}


#############################################################
## save_binStats
#' @description This function calculates prediction statistics for a set of binary prediction files and saves statistics to an output file.
#'
#' @param pred_files Character object including the file names containing binary predictions. If NULL, the function will use all the files ending in 'pred.txt' in the current directory (default NULL).
#' @param out_file Output file name.
#' @param digits Number of digits to round prediction statistics (default 3).
#'
#' @examples
#'
#'
#' save_binStats(inputFile = "binary_stats.txt", digits = 4)
#'
#' save_binStats(pred_files = c("Allantois_pred.txt", "Cardiomyocytes_pred.txt", "Endothelium_pred.txt")
#' , inputFile = "binary_stats.txt")
#'
#' @export
#'

save_binStats <- function(pred_dir = NULL, pred_files = NULL, out_file, digits = 3){
  if(is.null(pred_files)){
    message("Using files ending with pattern 'pred.txt'")
    
    if(is.null(pred_dir)){
      pred_files <- list.files(path = ".", pattern = "(.*)pred.txt$", full.names = TRUE)
    }else{
      pred_files <- list.files(path = pred_dir, pattern = "(.*)pred.txt$", full.names = TRUE)
    }
      
  }
  
  celltypes <- sub("_pred.txt", "", basename(pred_files))
  pred_li <- lapply(pred_files, read.table, header = TRUE, stringsAsFactors = FALSE)

  pred_stats <- lapply(1:length(pred_li)
                       , function(x) binStats_single(pred_li[[x]], celltypes[x]))
  pred_stats <- do.call("rbind", pred_stats)
  
  stat_cols <- vapply(pred_stats, is.numeric, FUN.VALUE = logical(1))
  pred_stats[, stat_cols] <- round(pred_stats[, stat_cols], digits = digits)
  
  write.table(x = pred_stats, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)
}
