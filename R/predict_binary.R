############################################################
#' add.missing.vars_xgb
#' @description 
#' Function to add columns with zero counts for motifs present in the model but not in the test set  
#' @param xgb.model XGBoost model
#' @param testSet Data frame containing a CRE test set 
#'
#' @export
add.missing.vars_xgb <- function(xgb.model, testSet){
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
#' extdata_path <- system.file("extdata", package = "BagOfMotifs")
#' motif_counts <- paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other_counts.txt")
#'
#'
#' predict_binary(motifs = motif_counts, xgb_model = paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other.rds"))
#' }
#' @export
predict_binary <- function(motifs, xgb_model, celltype_prefix = NULL
                           , pred = NULL, training_set = NULL, title = NULL)
{
  message("Reading classification model...")
  
  xgb <- readRDS(xgb_model)

  message("Saving XGBoost model in .bin format...")
  xgboost::xgb.save(xgb, gsub(".rds", ".bin", xgb_model))
  
  message(paste("Best tree:", xgb$best_iteration, "\n"))
  
  # Reading table of motif counts
  
  message("Reading motif counts matrix...")
  counts.tab <- read.table(file = motifs, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
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
  if (is.null(training_set)) {
    if(is.null(celltype_prefix)){
      stop("No celltype_prefix provided.")
    }
    training_set <- paste0(celltype_prefix, "_vs_Others_train.txt")
    message(paste("Saving training set to", training_set, "..."))
    write.table(x = motifs_train, file = training_set, quote = FALSE, sep = "\t")
  } else {
    message(paste("Saving training set to", training_set, "..."))
    write.table(x = motifs_train, file = training_set, quote = FALSE, sep = "\t")
  }
  
  motifs_test <- add.missing.vars_xgb(xgb, motifs_test)
  test_labels <- motifs_test$binary_celltype
  motifs_test <- motifs_test[, xgb$feature_names]  
  
  set.seed(123)
  y_pred <- predict(xgb, data.matrix(motifs_test), type = "response")
  predicted.class <- y_pred > 0.5 
  predicted.class <- gsub("TRUE", 1, predicted.class)
  predicted.class <- gsub("FALSE", 0, predicted.class)
  
  actual.vs.predicted <- data.frame(true_class = test_labels
                                    , predicted_class = predicted.class
                                    , prob = y_pred
                                    , stringsAsFactors = FALSE)
  
  rownames(actual.vs.predicted) <- rownames(motifs_test)
  
  if (is.null(pred)) {
    if(is.null(celltype_prefix)){
      stop("No celltype_prefix provided.")
    }
    pred <- paste0(celltype_prefix, "_vs_Others_pred.txt")
    message(paste("Saving predictions to", pred, "..."))
    write.table(x = actual.vs.predicted, file = pred, quote = FALSE, sep = "\t")
  } else {
    message(paste("Saving predictions to", pred, "..."))
    write.table(x = actual.vs.predicted, file = pred, quote = FALSE, sep = "\t")
  }
  
  message("Preparing ROC plots")
  require(cvAUC)
  require(pROC)
  require(ggplot2)
  
  
  roc_pred_tab <- roc(actual.vs.predicted$true_class, actual.vs.predicted$prob, direction = "<")
  roc_auc <- auc(roc_pred_tab)
  
  rocs.list <- list(roc_pred_tab)
  # Plot ROC curves
  p <- ggroc(rocs.list, size = 1) + theme(panel.border = element_blank(), panel.grid.major = element_blank()
                                          , panel.grid.minor = element_blank() #legend.position = "none" 
                                          , axis.line = element_line(colour = "black", size = 1)
                                          , legend.title = element_blank()
                                          , legend.key = element_blank()
                                          , legend.position = "none"
                                          , legend.key.width = unit(1.5, "cm")
                                          , panel.background = element_blank()
                                          , text = element_text(size = 23)
                                          , axis.text.x = element_text(colour = "black")
                                          , axis.text.y = element_text(colour = "black")
                                          , legend.text = element_text(size = 24)
                                          , axis.ticks = element_line(colour = "black")) +
    guides(linetype = guide_legend(override.aes = list(size = 3))) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", alpha = 0.8, color = "grey") + coord_equal() +
    scale_colour_manual(values = c("#B39EB5"), aesthetics = c("colour", "fill")) +
    labs(y = "Sensitivity", x = "Specificity") + 
    annotate("text", x = .5, y = .25, label = paste("AUC =", round(auc(roc_pred_tab), 2)))
  
  if (! is.null(title)){
	  p <- p + ggtitle(title)
  }
  
  return(p)
}

			   
			   
#############################################################
#' predict_binary_multi
#' @description Wrapper function for predict_binary. Will run 
#' on multiple files within specified directory
#'
#' @param inputMotif_dir  Directory of motif counts (named <celltype>_vs_Others.txt)
#' @param inputXGB_dir    Directory of XGBoost files (names <celltype>_vs_Others.rds)
#' @param outputTrain_dir Directory where to save training text output files
#' @param ncol				Number of columns of plots per page
#' @param nrow             Number of rows of plots per page
#' @param outputFile      Output pdf file name for plots 
#'
#' @examples 
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' motif_counts <- paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other_counts.txt")
#'
#'
#' predict_binary(motifs = motif_counts, xgb_model = paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other.rds")
#' }
#'
#' @export
#'
predict_binary_multi <- function(inputMotif_dir = NULL, inputXGB_dir = NULL, outputTrain_dir = NULL
                                 , pred_dir = NULL, ncol = 3, nrow = 3, 
                                 outputFile = NULL, width = 12, height = 12){
  if (is.null(outputFile )){
	  warning("Please provide output file name") 
	  return(-1)
  }
  require(cowplot)	
  
  fl.motifs <- list.files(path = inputMotif_dir, pattern = "*_vs_Others.txt")
  fl.xgb    <- list.files(path = inputXGB_dir, pattern = "*_vs_Others.rds")
  
  # Identify which entries have both a 
  candidates <- c(gsub(pattern = "\\.txt$", replacement = "", x = fl.motifs), gsub(pattern = "\\.rds$", replacement = "", x = fl.xgb))
  candidates <- candidates[which(duplicated(candidates))]
  
  allPlots <- list()
  
  for (i in 1: length(candidates)){
	  message(paste0( "Preparing ", candidates[i], " predict binary output") )
	  allPlots[[i]] <- suppressMessages(
		  predict_binary(motifs = paste0(inputMotif_dir,  "/", candidates[i], ".txt"),
				 xgb_model = paste0(inputXGB_dir,    "/", candidates[i], ".rds"), 
				 training_set  = paste0(outputTrain_dir, "/", candidates[i], "_train.txt"),
				 pred = paste0(pred_dir, "/", candidates[i], "_pred.txt"),
				 title = gsub(pattern = "_vs_Others", replacement = "", x = candidates[i])))
  }
	
  nsamples <- ncol * nrow
  currentSample <- 1
  pages <- ceiling(length(candidates) / nsamples)
  
  
  pdf(outputFile, width = width, height = height)
  for(i in 1:pages) {	
    print(plot_grid(plotlist=allPlots[currentSample:(nsamples*i)], ncol = ncol, nrow = nrow))
    currentSample <- currentSample + nsamples
  }
  dev.off()
  
  
  
}



#############################################################
#' prcurve_binary
#' @description Function to produce a precision/recall curve for binary prediction  
#' 
#' @param motifs  File containing motif counts
#' @param xgb_model File name of model
#' @param title Title for output plot (default NULL)
#' 
#' @examples 
#' \dontrun{
#' 
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' motif_counts <- paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other_counts.txt")
#' 
#' 
#' predict_binary(motifs = motif_counts, xgb_model = paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other.rds"))
#' }
#' @export
prcurve_binary <- function(motifs, xgb_model, title = NULL){
  message("Reading classification model...")
  
  xgb <- readRDS(xgb_model)
  # Reading table of motif counts
  
  message("Reading motif counts matrix...")
  counts.tab <- read.table(file = motifs, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
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
  
  motifs_test <- add.missing.vars_xgb(xgb, motifs_test)
  test_labels <- motifs_test$binary_celltype
  motifs_test <- motifs_test[, xgb$feature_names]  
  
  set.seed(123)
  y_pred <- predict(xgb, data.matrix(motifs_test), type = "response")
  predicted.class <- y_pred > 0.5 
  predicted.class <- gsub("TRUE", 1, predicted.class)
  predicted.class <- gsub("FALSE", 0, predicted.class)
  
  actual.vs.predicted <- data.frame(true_class = test_labels
                                    , predicted_class = predicted.class
                                    , prob = y_pred
                                    , stringsAsFactors = FALSE)
  
  rownames(actual.vs.predicted) <- rownames(motifs_test)

  message("Preparing PR curve")
  require(yardstick)
  require(cvAUC)
  require(pROC)
  require(ggplot2)
  
  actual.vs.predicted$true_class <- factor(actual.vs.predicted$true_class, levels = c(0,1))
  pr_curve.df <- yardstick::pr_curve(actual.vs.predicted, truth = true_class, prob
                                     , event_level = "second")
  
  pr_auc_val <- yardstick::pr_auc(actual.vs.predicted, truth = true_class, prob
                                  , event_level = "second", estimator = "binary")
  
  # plot precision recall curve
  p <- ggplot(pr_curve.df, aes(x = recall, y = precision)) +
    geom_path() +  theme(panel.border = element_blank(), panel.grid.major = element_blank()
                                         , panel.grid.minor = element_blank() #legend.position="none" 
                                         , axis.line = element_line(colour = "black", size = 1)
                                         , legend.title = element_blank()
                                         , legend.key = element_blank()
                                         , legend.position = "none"
                                         , legend.key.width = unit(1.5, "cm")
                                         , panel.background = element_blank()
                                         , text = element_text(size = 23)
                                         , axis.text.x = element_text(colour = "black")
                                         , axis.text.y = element_text(colour = "black")
                                         , legend.text = element_text(size = 24)
                                         , axis.ticks = element_line(colour = "black")) +
    guides(linetype = guide_legend(override.aes = list(size = 3))) +
    coord_equal() + scale_colour_manual(values = c("#B39EB5"), aesthetics = c("colour", "fill")) +
    labs(y = "Precision", x = "Recall") + 
    annotate("text", x = .5, y = .25, label = paste("AUC =", round(pr_auc_val$.estimate, 2)))
  
  
  if (! is.null(title)){
	  p <- p + ggtitle(title)
  }
  
  return(p)
}



#############################################################
#' prcurve_multi
#' @description Wrapper function for prcurve_binary. Will run 
#' on multiple files within specified directory
#'
#' @param inputMotif_dir  Directory of motif counts (named <celltype>_vs_Others.txt)
#' @param inputXGB_dir    Directory of XGBoost files (names <celltype>_vs_Others.rds)
#' @param ncol				Number of columns of plots per page
#' @param nrow             Number of rows of plots per page
#' @param outputFile      Output pdf file name for plots 
#' @param width       Width of the graphics region in pdf output in inches (default 12) 
#' @param height      Height of the graphics region in pdf output in inches (default 12)
#'
#' @examples 
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' motif_counts <- paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other_counts.txt")
#'
#'
#' predict_binary(motifs = motif_counts, xgb_model = paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other.rds")
#' }
#'
#' @export
#'
prcurve_multi <- function(inputMotif_dir = NULL, inputXGB_dir = NULL
                                 , ncol = 3, nrow = 3, 
                                 outputFile = NULL, width = 12, height = 12){
  if (is.null(outputFile )){
	  warning("Please provide output file name")
	  return(-1)
  }
  require(cowplot)	
  
  fl.motifs <- list.files(path = inputMotif_dir, pattern = "*_vs_Others.txt")
  fl.xgb <- list.files(path = inputXGB_dir, pattern = "*_vs_Others.rds")
  
  # Identify which entries have both a 
  candidates <- c(gsub(pattern = "\\.txt$", replacement = "", x = fl.motifs), gsub(pattern = "\\.rds$", replacement = "", x = fl.xgb))
  candidates <- candidates[which(duplicated(candidates))]
  
  allPlots <- list()
  
  for (i in 1: length(candidates)){
	  message(paste0( "Preparing ", candidates[i], " PR curve"))
	  allPlots[[i]] <- suppressMessages(
		  prcurve_binary(motifs = paste0(inputMotif_dir,  "/", candidates[i], ".txt"),
				 xgb_model = paste0(inputXGB_dir,    "/", candidates[i], ".rds"), 
				 title = gsub(pattern = "_vs_Others", replacement = "", x = candidates[i])
      ))
  }
  nsamples <- ncol * nrow
  currentSample <- 1
  pages <- ceiling(length(candidates) / nsamples)
  
  
  pdf(outputFile, width = width, height = height)
  for(i in 1:pages){
    print(plot_grid(plotlist = allPlots[currentSample:(nsamples * i)], ncol = ncol, nrow = nrow))
    currentSample <- currentSample + nsamples
  }
  dev.off()
  
}

