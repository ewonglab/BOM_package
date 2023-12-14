#############################################################
#' binModel
#' @description Function to produce a matrix of motif frequency. Motif instances will be filtered using the q-value threshold provided. 
#' The regions annotated to the target cell type/condition will be coded as 1. Background regions will be coded as 0.
#' The outpud matrix will contain a balanced number of positive (target) and negative (background) instances.
#' 
#' @param target_ct  Name of the target cell type/condition. If target_ct is NULL, the function will produce a matrix for each cell type in the input directory. 
#' @param data_path Path to the input data directory.
#' @param qval_thresh q-value threshold for motif filtering. Default to 0.5 (q-value <= 0.5).
#' @param outDir Name of directory to save output files. Output files will be named cellType_vsOthers.
#' @param nthreads The number of threads to use
#' 
#' @examples 
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' data_path <- paste0(extdata_path, "/tutorial/motifs") 
#' 
#' binModel(data_path=data_path, qval_thresh=0.5, outDir='results/', target_ct= 'Cardiomyocytes')
#' }
#' @export
binModel <- function(data_path, qval_thresh, outDir, target_ct=NULL ,nthreads=1)
{
	require(foreach)
  # Set up multiple workers
  system.name <- Sys.info()['sysname']
  new_cl <- FALSE
  if (system.name == "Windows") {
    new_cl <- TRUE
    cluster <- parallel::makePSOCKcluster(rep("localhost", nthreads))
    doParallel::registerDoParallel(cluster)
  } else {
    doParallel::registerDoParallel(cores=nthreads)
  }





	# Check requested threads is available
#	coresAvailable <- parallel::detectCores()
#	if (coresAvailable < nthreads)
#	{	nthreads <- coresAvailable }
#	cl <- parallel::makeCluster(nthreads)


	# list directories containing FIMO output
	directories <- list.dirs(path = data_path, full.names = T, recursive=F)
	celltypes <- basename(directories)
  
	counts_files <- list.files(path=data_path, full.names =T, recursive=T, pattern='fimo.tsv*')
	# Read results of motif search
	message(paste0("Reading input data from ",data_path,".\nThere are ", length(directories), " directories found."))
	
	read_and_update <- function(fn, pb) 
		{  
			cat("=")
			as.data.frame(read.table(fn, sep='\t', header=TRUE))
		}

	counts <- foreach::foreach(thisSample=counts_files) %dopar% {
		read_and_update(thisSample)
	}
#	suppressWarnings({
#		counts <- parallel::parLapply(cl, counts_files, read_and_update)
#	})
#	close(pb)
  
	counts <- lapply(counts, as.data.frame)
  for(i in 1:length(counts))
  {
    if (nrow(counts[[i]]) == 0)
    {
      stop(paste0("No data for ",countfiles[i]))
    }
  }


	counts <- lapply(counts, function(x) x[x[,9] <= qval_thresh,])
	names(counts) <- celltypes

  
	## count the number of unique CREs per condition
	n_CREs_by_ct <- unlist(lapply(counts, function(x) length(unique(x$sequence_name))))
	names(n_CREs_by_ct) <- celltypes
	message("Writing output")
	if (! is.null(target_ct)) # Do one requested comparison
	{   message(paste0("Doing selected comparison on ", target_ct))
		## define the number of CREs from target and background conditions
		if(! target_ct %in% names(n_CREs_by_ct)){
			stop(paste(target_ct, "is not among the contexts provided. Please check the spelling and case."))
		}
		binModel_oneVsOthers(target_ct, counts, n_CREs_by_ct, celltypes, outDir)
	}
	else  # Do all comparisons
	{	message("Processing all cell types")

	#	parallel::parLapply(cl, celltypes, binModel_oneVsOthers, counts, n_CREs_by_ct, celltypes, outDir)
#browser()
	outputforDebug <-foreach::foreach(thisCellType=celltypes) %dopar% {
		binModel_oneVsOthers(thisCellType, counts, n_CREs_by_ct, celltypes, outDir)
	} 

		# Train
		
		
		message("training....")
		inputData <- paste0(outDir,"/",celltypes,"_vs_Others.txt")
		outputFileNames <- paste0(outDir,"/",celltypes,"_vs_Others.rds")
		
	#	mapply(train_binary, inputData, save_name = outputFileNames, nthread = nthreads, verbose = 0)
		for ( i in 1:length(inputData))
		{
			message(paste0("Preparing training for ",celltypes[i]))
			train_binary(input_data=inputData[i], save_name=outputFileNames[i], 
								early_stopping_rounds=100,   verbose=0, nthread=nthreads)

		}
		
		
		
	}
#	parallel::stopCluster(cl)
	
	
	if (new_cl) { ## Shut down cluster if on Windows
    ## stop cluster
		parallel::stopCluster(cluster)
	}
	
	
	
	message("Complete")
}

################################################################################
## binModel_oneVsOthers
#'
#' @param target_ct  Target cell type name
#' @param n_CREs_by_ct  number of cis regulatory elements by cell type
#' @param celltypes   number of cell types
#' @param outDir out directory
#'
#'
binModel_oneVsOthers <- function(target_ct, counts, n_CREs_by_ct, celltypes, outDir)
{
  n <- n_CREs_by_ct[target_ct]
  n_bkg <- round(n/(length(celltypes) -1))
  
  # check whether the number of CREs of every context is enough to make a balanced set

  if(n_bkg > min(n_CREs_by_ct)){
    warning("The number of CREs for a background set is lower than the required number. Reducing the number of CREs from each context.")
    n_bkg <- min(n_CREs_by_ct)
    n <- n_bkg * (length(celltypes) -1)
    tmp <- counts[[target_ct]]
    set.seed(123)
    pos_sample <- sample(x = unique(tmp$sequence_name), size = n, replace = F)
    counts[[target_ct]] <- tmp[tmp$sequence_name %in% pos_sample,]
  }
  
  ###      Preparing negative (background) set
  
  negative_set <- lapply(celltypes[celltypes != target_ct], function(ct) {
    tmp <- counts[[ct]]

    set.seed(123)
    ct_sample <- sample(unique(tmp$sequence_name), size = n_bkg, replace = FALSE)
    tmp[tmp$sequence_name %in% ct_sample, ]
  })
  
  names(negative_set) <- celltypes[celltypes!=target_ct]
  
  for (i in seq_along(negative_set)) {
    negative_set[[i]]$celltype <- celltypes[celltypes != target_ct][i]
  }
  
  negative_set.df <- do.call("rbind", negative_set)
  
  negative_set.df <- negative_set.df[,c("motif_id", "sequence_name", "celltype")]
  tmp <- as.data.frame(table(negative_set.df[,1:2]))
  
  tmp <- tidyr::spread(tmp, motif_id, Freq)
  tmp <- merge(tmp, unique(negative_set.df[,2:3]), by="sequence_name")
#browser()
  if (length(unique(tmp$sequence_name)) != nrow(tmp))
  { 	warning("Unexpectedly Input data has duplicate peaks for ",target_ct,"! (when generating background)") 
		rownames(tmp) <- make.unique(tmp$sequence_name)
  }
  else
  {
	rownames(tmp) <- tmp$sequence_name 
  } 
  tmp$sequence_name <- NULL

  ### Preparing positive (target) set
  
  positive <- counts[[target_ct]]
  positive$celltype <- target_ct
  positive <- positive[,c("motif_id", "sequence_name", "celltype")]
  tmp2 <- as.data.frame(table(positive[,1:2]))
  tmp2 <- tidyr::spread(tmp2, motif_id, Freq)
  tmp2 <- merge(tmp2, unique(positive[,2:3]), by="sequence_name")
  #enhancer IDs as rownames
  if (length(unique(tmp2$sequence_name)) != nrow(tmp2))
  { 	warning("Unexpectedly Input data has duplicate peaks for ",target_ct,"! (when generating positive)") 
		rownames(tmp2) <- make.unique(tmp2$sequence_name)
  }
  else
  {
	rownames(tmp2) <- tmp2$sequence_name 
  } 
  
  ## Combine target and background sets

  final_set <- dplyr::bind_rows(tmp, tmp2)
  final_set[is.na(final_set)] <- 0
  final_set <- final_set[,c(colnames(final_set)[colnames(final_set)!="celltype"],"celltype")]
  final_set$binary_celltype <- ifelse(final_set$celltype==target_ct, 1, 0)
  
#  if()
  message(paste("Saving matrix of motif counts for binary classification...\n"))
  write.table(x = final_set, file = paste0(outDir,"/",target_ct,"_vs_Others.txt"), quote = F, sep = '\t')
  
  message("Content of output table:\n")
  out_content <- as.data.frame(table(final_set$binary_celltype))
  colnames(out_content) <- c("Context", "Number of elements")
  print(out_content)
  message(paste0("Processed ", target_ct))
}

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
#' @param training Proportion of data that will be split for training and test. 
#' Default is 0.6 which leaves 0.2 and 0.2 for validation and testing respectively.
#'      
#' @examples 
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' motif_counts <- paste0(extdata_path, "/tutorial/motifs/Cardiomyocytes_vs_other_counts.txt")
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
                         , callbacks = list()
						 , training = 0.6)
{
	if ((training > 1) | (training < 0) )
	{ error("Parameter training has to be between 0 and 1") }

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
  if (verbose)
  {  cat("Reading motif counts...\n")  }
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
  if (verbose)
  {  message("Splitting data into training, validation and test sets...\n")  }
  set.seed(123)
  motifs_split <- rsample::initial_split(counts.tab, prop = training)
  motifs_train <- rsample::training(motifs_split)
  motifs_test <- rsample::testing(motifs_split)
  
  set.seed(123)
  motifs_split2 <- rsample::initial_split(motifs_test, prop = (1-training)/2)
  motifs_val <- rsample::training(motifs_split2)
  motifs_test <- rsample::testing(motifs_split2)

  # Removing non-variable motifs from training set
  motifs_train.sd <- apply(motifs_train, 2, sd)
  motifs_train <- motifs_train[, names(which(motifs_train.sd != 0)), drop = F]
  motifs_val <- motifs_val[,colnames(motifs_train), drop = F]
  
  # Prepare training and validation DMatrix objects
  dtrain <- xgboost::xgb.DMatrix(label = as.numeric(motifs_train$binary_celltype)
                        , data = as.matrix(motifs_train[, colnames(motifs_train)[colnames(motifs_train)!="binary_celltype"]]))
  dvalid <- xgboost::xgb.DMatrix(label = as.numeric(motifs_val$binary_celltype)
                    , data = as.matrix(motifs_val[, colnames(motifs_val)[colnames(motifs_val)!="binary_celltype"]]))
  
  # params$data <- dtrain
  watchlist <- list(train = dtrain, validation = dvalid) 
  
  if (verbose)
  { message("Training binary classification model...\n") }
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
	if(verbose)
	{	message("Saving model...\n")  }
  
	#xgb.save(model, params$save_name)
	saveRDS(model, save_name)
}

