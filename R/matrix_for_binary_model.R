#############################################################
#' matrix_binModel
#' @description Function to produce a matrix of motif frequency. Motif instances will be filtered using the q-value threshold provided. 
#' The regions annotated to the target cell type/condition will be coded as 1. Background regions will be coded as 0.
#' The outpud matrix will contain a balanced number of positive (target) and negative (background) instances.
#' 
#' @param target_ct  Name of the target cell type/condition
#' @param data_path Path to the data directory
#' @param qval_thresh q-value threshold for motif filtering
#' @param outDir Name of directory to save output files. Output files will be named cellType_vsOthers.
#' @param nthreads The number of threads to use
#' 
#' @examples 
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' data_path <- paste0(extdata_path, "/tutorial/motifs") 
#' binPredictions <- paste0(extdata_path, "/tutorial/Cardiomyocytes_vs_other_pred.txt")
#' 
#' mat <- matrix_binModel(data_path=data_path, qval_thresh=0.5, outDir='results/', target_ct= 'Cardiomyocytes')
#' }
#' @export
matrix_binModel <- function(data_path, qval_thresh, outDir, target_ct=NULL ,nthreads=1)
{
	# Check requested threads is available
	coresAvailable <- parallel::detectCores()
	if (coresAvailable < nthreads)
	{	nthreads <- coresAvailable }
	cl <- parallel::makeCluster(nthreads)


	# list directories containing FIMO output
	directories <- list.dirs(path = data_path, full.names = T, recursive=F)
	celltypes <- basename(directories)
  
	counts_files <- list.files(path=data_path, full.names =T, recursive=T, pattern='fimo.tsv*')
	# Read results of motif search
	message(paste0("Reading input data from ",data_path,".\nThere are ", length(directories), " directories found."))
	
	read_and_update <- function(fn, pb) 
		{  
			setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
			cat("=")
			as.data.frame(read.table(fn, sep='\t', header=TRUE))
		}
	
	pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(celltypes), # Maximum value of the progress bar
                     style = 3,    				# Progress bar style (also available style = 1 and style = 2)
                     width = 50,   			# Progress bar width. Defaults to getOption("width")
                     char = "=")   			# Character used to create the bar	
	
	suppressWarnings({
    #counts <- lapply(counts_files, data.table::fread)
		counts <- parallel::parLapply(cl, counts_files, read_and_update, pb)
		#lapply(counts_files, read.table,  sep="\t")
	})
	close(pb)
  
	counts <- lapply(counts, as.data.frame)
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
	{	message("Processign all cell types")
#		pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
#                     max = length(celltypes), # Maximum value of the progress bar
#                     style = 3,    				# Progress bar style (also available style = 1 and style = 2)
#                     width = 50,   			# Progress bar width. Defaults to getOption("width")
#                     char = "=")   			# Character used to create the bar
#		for ( i in 1:length(celltypes))
#		{	
#			binModel_oneVsOthers(celltypes[i], n_CREs_by_ct, celltypes, outDir)
#			setTxtProgressBar(pb, i)
#		}
#		close(pb)
		parallel::parLapply(cl, celltypes, binModel_oneVsOthers, counts, n_CREs_by_ct, celltypes, outDir)
		
	}
	parallel::stopCluster(cl)
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
  rownames(tmp) <- tmp$sequence_name
  tmp$sequence_name <- NULL

  ### Preparing positive (target) set
  
  positive <- counts[[target_ct]]
  positive$celltype <- target_ct
  positive <- positive[,c("motif_id", "sequence_name", "celltype")]
  tmp2 <- as.data.frame(table(positive[,1:2]))
  tmp2 <- tidyr::spread(tmp2, motif_id, Freq)
  tmp2 <- merge(tmp2, unique(positive[,2:3]), by="sequence_name")
  #enhancer IDs as rownames
  rownames(tmp2) <- tmp2$sequence_name
  tmp2$sequence_name <- NULL
  
  ## Combine target and background sets

  final_set <- dplyr::bind_rows(tmp, tmp2)
  final_set[is.na(final_set)] <- 0
  final_set <- final_set[,c(colnames(final_set)[colnames(final_set)!="celltype"],"celltype")]
  final_set$binary_celltype <- ifelse(final_set$celltype==target_ct, 1, 0)
  
#  if()
  message(paste("Saving matrix of motif counts for binary classification...\n"))
  write.table(x = final_set, file = paste0(outDir,"/",target_ct,"_vs_Others"), quote = F, sep = '\t')
  
  message("Content of output table:\n")
  out_content <- as.data.frame(table(final_set$binary_celltype))
  colnames(out_content) <- c("Context", "Number of elements")
  print(out_content)
  message(paste0("Processed ", target_ct))
}

