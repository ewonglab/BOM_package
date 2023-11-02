#############################################################
#' train_binary
#' @description Function to train a binary classification model. 
#' 
#' @param input_data  File containing the input matrix of motif counts
#' @param nrounds Number of boosting rounds (default: 10000)
#' @param eta Learning rate (default: 0.01)
#' @param training_set File name to output training counts (optional)
#' @param max_depth Maximum tree depth (default: 6)
#' @param subsample Subsample ratio of the training instances (default: 0.5)
#' @param colsample_bytree Subsample ratio of columns when constructing each tree (default: 0.5)
#' @param objective Objective function (default: binary:logistic)
#' @param early_stopping_rounds Perform early stopping if no improvement for this many rounds (default: NULL)
#' @param nthread Number of parallel threads (default: 1)
#' @param eval_metric Evaluation metric (default: error)
#' @param maximize Whether to maximize the evaluation metric (default: FALSE)
#' @param save_period Save the model for every given number of rounds (default: NULL)
#' @param save_name Name of the saved model file (default: xgboost.model)
#' @param feval Customized evaluation metric (default: NULL)
#' @param Verbose How much details on the progress to print (default: 1)
#' @param print_every_n Print evaluation messages each n-th iterations (default: 1)
#'      
#' @examples 
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' motif_counts <- paste0(extdata_path, "/tutorial/Cardiomyocytes_vs_other_counts.txt")
#' 
#' 
#' train_binary(input_data = motif_counts, save_name="Cardiomyocytes_vs_other.rds"
#' , early_stopping_rounds = 100, print_every_n = 100, nthread = 4)
#' }
#' @export
train_binary <- function(input_data = NULL, nrounds = 10000
                         , eta = 0.01, max_depth = 6, subsample = 0.5
                         , colsample_bytree = 0.5, objective = "binary:logistic"
                         , early_stopping_rounds = NULL
                         , nthread = 1, eval_metric = "error", maximize = F
                         , params = list(), feval = NULL, verbose = 1
                         , print_every_n = 1L, save_period = NULL
                         , save_name = "xgboost.model", xgb_model = NULL
                         , callbacks = list())
{
  nrounds <- as.integer(nrounds)
  max_depth <- as.integer(max_depth)
  
  if(!is.null(early_stopping_rounds)){
    early_stopping_rounds <- as.integer(early_stopping_rounds)
  }
  
  nthread <- as.integer(nthread)
  
  if(!is.null(save_period)){
    save_period <- as.integer(save_period)
  }
  
  if(!is.null(nthread)){
    nthread <- as.integer(nthread)
  }
  
  eta <- as.numeric(eta)
  subsample <- as.numeric(subsample)
  colsample_bytree <- as.numeric(colsample_bytree)
  maximize <- as.logical(maximize)
  
  # Reading table of motif counts
  cat("Reading motif counts...\n")
  counts.tab <- read.table(file = input_data, header =T
                           , stringsAsFactors = F, sep = '\t')
  counts.tab$celltype <- NULL
  
  counts.tab.NAs <- sapply(counts.tab, function(x) sum(is.na(x)))
  if(any(counts.tab.NAs) > 0){
    warning("NAs present in matrix of motif counts...\n")  
  }
  
  # Binary label as numeric
  counts.tab$binary_celltype <- as.numeric(counts.tab$binary_celltype)
  
  # Split dataset into training, validation and test sets
  message("Splitting data into training, validation and test sets...\n")
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
  motifs_train <- motifs_train[, names(which(motifs_train.sd != 0)), drop = F]
  motifs_val <- motifs_val[,colnames(motifs_train), drop = F]
  
  # Prepare training and validation DMatrix objects
  dtrain <- xgb.DMatrix(label = as.numeric(motifs_train$binary_celltype)
                        , data = as.matrix(motifs_train[, colnames(motifs_train)[colnames(motifs_train)!="binary_celltype"]]))
  dvalid <- xgb.DMatrix(label = as.numeric(motifs_val$binary_celltype)
                        , data = as.matrix(motifs_val[, colnames(motifs_val)[colnames(motifs_val)!="binary_celltype"]]))
  
  # params$data <- dtrain
  watchlist <- list(train = dtrain, validation = dvalid) 
  
  
  message("Training binary classification model...\n")
  set.seed(123) 
  model <- xgboost::xgb.train(
    data = dtrain,
    nrounds = nrounds,
    watchlist = watchlist,
    objective = objective,
    eta = eta,
    max_depth = max_depth,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    nthread = nthread,
    eval_metric = eval_metric,
    params = params,
    feval = feval,
    verbose = verbose,
    print_every_n = print_every_n,
    early_stopping_rounds = early_stopping_rounds,
    maximize = maximize,
    save_period = save_period,
    save_name = save_name,
    xgb_model = xgb_model,
    callbacks = callbacks
  )
  
  message("Saving model...\n")
  
  #xgb.save(model, params$save_name)
  saveRDS(model, save_name)
}

