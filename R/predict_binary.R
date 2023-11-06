############################################################
#' add.missing.vars_xgb
#' @description 
#' Function to add columns with zero counts for motifs present in the model but not in the test set  
#' @param xgb.model XGBoost model
#' @param testSet Data frame containing a CRE test set 
#'
#' @export
add.missing.vars_xgb <- function(xgb.model, testSet)
  {
  missing.vars <- setdiff(xgb.model$feature_names, colnames(testSet))
  testSet[, missing.vars] <- lapply(missing.vars, function(var) as.integer(0))
  return(testSet)
  }


#############################################################
#' predict_binary
#' @description Function to predict a set of test CREs using a binary classification model  
#' 
#' @param motifs  File containing motif counts
#' @param xgb_model File name of model
#' @param pred File name to output predicted values (default "predictions.txt")
#' @param training_set File name to output training counts (optional)
#' 
#' @examples 
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' motif_counts <- paste0(extdata_path, "/tutorial/Cardiomyocytes_vs_other_counts.txt")
#' 
#' 
#' predict_binary(motifs = motif_counts, xgb_model = paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other.rds")
#' }
#' @export
predict_binary <- function(motifs, xgb_model, training_set = NULL, pred = "predictions.txt")
  {
  message("Reading classification model...")
  
  xgb <- readRDS(xgb_model)
  
  message("Saving XGBoost model in .bin format...")
  xgboost::xgb.save(xgb, gsub(".rds", ".bin", xgb_model))
  
  message(paste("Best tree:", xgb$best_iteration, "\n"))
  
  # Reading table of motif counts
  
  message("Reading motif counts matrix...")
  counts.tab <- read.table(file = motifs, header =T, stringsAsFactors = F, sep = '\t')
  counts.tab$celltype <- NULL
  
  counts.tab.NAs <- sapply(counts.tab, function(x) sum(is.na(x)))
  
  if(any(counts.tab.NAs) > 0){
    warning("NAs present in input matrix...")  
  }
  
  # Binary label as numeric
  counts.tab$binary_celltype <- as.numeric(counts.tab$binary_celltype)
  
  # Split dataset into training, validation and test sets
  set.seed(123)
  motifs_split <- rsample::initial_split(counts.tab, prop = .6)
  motifs_train <- rsample::training(motifs_split)
  motifs_test <- rsample::testing(motifs_split)
  
  set.seed(123)
  motifs_split2 <- rsample::initial_split(motifs_test, prop = .5)
  motifs_val <- rsample::training(motifs_split2)
  motifs_test <- rsample::testing(motifs_split2)
  
  # Removing non-variable motifs from training set
  motifs_train.sd <- apply(motifs_train, 2, sd)
  motifs_train <- motifs_train[, names(which(motifs_train.sd != 0))]
  
  # Save training set if a file name is provided
  if (!is.null(training_set)) {
    message(paste("Saving training set to", training_set, "..."))
    write.table(x = motifs_train, file = training_set, quote = FALSE, sep ='\t')
  }
  
  motifs_test <- add.missing.vars_xgb(xgb, motifs_test)
  test_labels <- motifs_test$binary_celltype
  motifs_test <- motifs_test[,xgb$feature_names]  
  
 set.seed(123)
  y_pred <- predict(xgb, data.matrix(motifs_test), type="response")
  predicted.class <- y_pred > 0.5 
  predicted.class <- gsub("TRUE", 1, predicted.class)
  predicted.class <- gsub("FALSE", 0, predicted.class)
  
  actual.vs.predicted <- data.frame(true_class = test_labels
                                    , predicted_class = predicted.class
                                    , prob = y_pred
                                    , stringsAsFactors = F)
  
  rownames(actual.vs.predicted) <- rownames(motifs_test)
  
  message(paste("Saving predicted values to", pred, "..."))
  write.table(x = actual.vs.predicted, file = pred, quote = F)

}
