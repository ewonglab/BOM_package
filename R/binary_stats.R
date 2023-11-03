############################################################
### MCC
#' Function to calculate Matthew's correlation coefficient (MCC)  
#' @param TP Number of True Positive instances
#' @param TN Number of True Negative instances
#' @param FP Number of False Positive instances
#' @param FN Number of False Negative instances
#'
MCC <- function(TP, TN, FP, FN){
  TP <- as.numeric(TP)
  TN <- as.numeric(TN)
  FP <- as.numeric(FP)
  FN <- as.numeric(FN)
  mcc_val <- ((TP * TN) - (FP * FN)) / 
    sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  return(mcc_val)
}

#############################################################
## binaryStats
#'
#' @param inputFile  input binary model predictions filename
#'
#' @examples 
#'
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' binPredictions <- paste0(extdata_path, "/tutorial/Cardiomyocytes_vs_other_pred.txt")
#' 
#' 
#' binaryStats(inputFile = binPredictions)
#
#' @export
binaryStats <- function(inputFile)
{
  message("Reading predicted values...\n")
  
  pred_tab <- read.table(file = inputFile, header = T, stringsAsFactors = F)
  TP <- nrow(pred_tab[pred_tab$true_class == 1 & pred_tab$predicted_class == 1, ])
  TN <- nrow(pred_tab[pred_tab$true_class == 0 & pred_tab$predicted_class == 0, ])
  FP <- nrow(pred_tab[pred_tab$true_class == 0 & pred_tab$predicted_class == 1, ])
  FN <- nrow(pred_tab[pred_tab$true_class == 1 & pred_tab$predicted_class == 0, ])

  recall <- TP/(TP + FN)
  precision <- TP/(TP + FP)
  f1 <- 2*((precision*recall)/(precision+recall))
  pred_tab$true_class <- factor(pred_tab$true_class, levels = c(0,1))
  
  pr_auc_val <- yardstick::pr_auc(pred_tab, truth = true_class, prob
                                  , event_level = "second", estimator = "binary")
  mcc <- MCC(TP, TN, FP, FN)
  Accuracy <- nrow(pred_tab[pred_tab$predicted_class==pred_tab$true_class,])/nrow(pred_tab)
  mcc <- MCC(TP, TN, FP, FN)
  
  
  auROC <- round(cvAUC::AUC(pred_tab$prob, pred_tab$true_class),4)
  pr_auc_val <- round(pr_auc_val$.estimate,4)
  Accuracy <- round(Accuracy, 4)
  f1 <- round(f1, 4)
  recall <- round(recall, 4)
  precision <- round(precision, 4)
  mcc <- round(mcc, 4)
  
  ## Print binary prediction statistics
  cat(paste("auROC:", auROC, '\n'))
  cat(paste("auPR:", pr_auc_val, '\n'))
  cat(paste("Accuracy:", Accuracy, '\n'))
  cat(paste("F1 score:", f1, '\n'))
  cat(paste("Recall:", recall, '\n'))
  cat(paste("Precision:", precision, '\n'))
  cat(paste("MCC:", mcc, '\n'))
  
}
