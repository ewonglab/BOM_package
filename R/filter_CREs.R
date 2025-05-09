############################################################
#' adjust_CREs
#' Function to trim (or extend) CREs 
#' @param x data frame of genomic coordinates. BED file format. 
#' @param N Number of base pairs - used to extend coordinates stored in x.
#' @param chom_sizes  dataframe of chromosome sizes
#'
adjust_CREs <- function(x, N, chrom_sizes){
  
  x$width <- with(x, V3 - V2)
  x$centre <- with(x, V2 + round(width/2))
  x$start <- round(x$centre - N/2)
  x$end <- round(x$centre + N/2)

  # Add chromosome size
  x$ord <- 1:nrow(x)
  x <- merge(x, chrom_sizes, by.x = "V1", by.y = "chr")

  # Look for CREs at the extremes of chromosomes
  if(any(x$start < 0) | any(x$end > x$chr_size)){
    warning(paste("Warning: some CREs are shorter than ", paste0(N, "bp")
                  , "because they are at the extremes of chromosomes"))
    if(any(x$start < 0)){
      # special case 1: negative start position 
      x$start <- ifelse(x$start < 0, 0, x$start)
    }
    if(any(x$end > x$chr_size)){
      # special case 2: CRE end position exceeds chromosome length
      x$end <- with(x, ifelse(end > chr_size, chr_size, end))
    }
  }
  
  # Return to original order
  x <- x[with(x, order(ord)), ]
  x <- x[, colnames(x)[!colnames(x) %in% c("chr_size", "centre", "width", "V2", "V3", "ord")]]
  x <- x[, c("V1", "start", "end", setdiff(colnames(x), c("V1", "start", "end")))]
  colnames(x)[2:3] <- c("V2", "V3")
  return(x)
  
}




#############################################################
#' filterCREs
#' @description Function to filter a bed file of Cis 
#' Regulatory Elements (CRE). Will create output bed files.
#'
#' @param inputBedFile  input bed filename
#' @param u number of basepairs upstream of transcriptional start sites
#' @param d number of basepairs downstream of transcriptional start site
#' @param nbp  number of base pairs that each CRE will be centered and adjusted to. i.e. each CRE will be adjusted to this width. 
#' @param chrSizesFile  File name of chromsome sizes ( tab 
#' delimited file with no header: first column chromosome ID, second column size)
#' @param keep_proximal  boolean. Only keep proximal regions (defined by u and d) relative to transcription start sites
#' @param remove_proximal  Boolean: Remove proximal regions (defined by u and d) relative to transcriptional start sites
#' @param non_exonic Boolean - defines if non exonic regions should ONLY be considered.
#' @param out_bed   filename of output bed file after filtering.
#' @param inputBedZeroBased Boolean (default TRUE). Add 1 to start position of input bed file.
#' @param celloutputDir  directory name (will create) where individual BED files will be populated. If NULL this step is ommited.
#' @param minCellPercent Each cell population must make up this percent of all cells (default 1 percent). Will notify what cells are removed
#' @param ovr_dir Boolean (default FALSE). Whether to overwrite output directory
#'
#' @examples 
#' \dontrun{
#' extdata_path <- system.file("extdata", package = "BagOfMotifs")
#' input.bed <- paste0(extdata_path, "/tutorial/mouseE8.25_peaks.bed")
#' annot.file <- paste0(extdata_path, "/Mus_musculus.GRCm38.92.gtf.gz")
#' mouseChrSizes <- paste0(extdata_path, "/mouse_sizes_primary_genome_tab.txt")
#'
#'
#' filterCRE(inputBedFile = input.bed, annotFile = annot.file, 
#'                     remove_proximal = TRUE,  non_exonic = TRUE,
#' 					   u =  1000, d=1000, nbp=500,
#'                     chrSizesFile = mouseChrSizes, 
#'                     out_bed =  'mouseE8.25_peaks_filt.bed')
#' }
#' @export
filterCREs <- function(inputBedFile = NULL,
                       annotFile = NULL,
                       chrSizesFile = NULL,
                       u = NULL,
                       d = NULL,
                       nbp = 500,
                       keep_proximal = FALSE,
                       remove_proximal = FALSE,
                       non_exonic = FALSE,
                       out_bed = NULL, inputBedZeroBased = TRUE,
                       celloutputDir = NULL, minCellPercent = 1,
                       ovr_dir = FALSE){
  if (is.null(celloutputDir)){
    stop("celloutputDir is not set.")
  } else if (file.exists(celloutputDir) & !ovr_dir){
      stop(paste("The directory", celloutputDir, "already exist. Set ovr_dir to TRUE"))
  } else {
      suppressWarnings(dir.create(celloutputDir))
  }
  
  message("Reading CREs...\n")
  addToBed <- 0
  if (inputBedZeroBased == TRUE){
	  addToBed <- 1
  } 
  
  cres <- read.table(file = inputBedFile, header = FALSE, stringsAsFactors = FALSE, sep = '\t')
  
  # Clean up 4th column by removing spaces an slashes
  cres$V4 <- sub("/", "_", cres$V4)
  cres$V4 <- sub(" ", "_", cres$V4)
  
  cres_gr <- with(cres, GenomicRanges::GRanges(V1, IRanges::IRanges(V2 + addToBed, V3), cellType = cres$V4))
  
  if(keep_proximal | remove_proximal | non_exonic){
    
    message("Reading genome annotation...\n")
    txdb_obj <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(file = annotFile, format = "gtf"))
    
    matchingChromosomes <- intersect(GenomeInfoDb::seqlevels(txdb_obj), GenomeInfoDb::seqlevels(cres_gr))
    if (length(matchingChromosomes) == 0){
      # No chromosomes align - Change style and try again
      GenomeInfoDb::seqlevelsStyle(cres_gr) <- "NCBI"
      matchingChromosomes <- intersect(GenomeInfoDb::seqlevels(txdb_obj), GenomeInfoDb::seqlevels(cres_gr))
      if (length(matchingChromosomes) == 0){
        warning("Seqence levels do not match to reference")
        warning(seqlevels(txdb_obj))
        warning(seqlevels(cres_gr))
      }
    }
    
    
    if(non_exonic){
      message("Removing exonic regions...\n")
      exons <- GenomicFeatures::exons(txdb_obj)
      x <- as.data.frame(GenomicRanges::findOverlaps(cres_gr, exons))
      if(nrow(x) > 0){
			cres_gr <- cres_gr[-unique(x$queryHits)]
      } 
    }
    
    if(keep_proximal & remove_proximal){
      stop(paste("'keep_proximal' and 'remove_proximal' are mutually exclusive"))
    }
    if(keep_proximal){
      message("Keeping proximal regions to TSSs...\n")
      proximal <- GenomicFeatures::promoters(x = txdb_obj, upstream = u, downstream = d)
      x <- as.data.frame(GenomicRanges::findOverlaps(cres_gr, proximal))
      if(nrow(x) > 0){
        cres_gr <- cres_gr[unique(x$queryHits)]
      } else {
        stop(paste("No remaining CREs"))
      }
    }
	  
    if(remove_proximal){
      message("Removing proximal regions to TSSs...\n")
      proximal <- GenomicFeatures::promoters(x = txdb_obj, upstream = u, downstream = d)
      x <- as.data.frame(GenomicRanges::findOverlaps(cres_gr, proximal))
      if(nrow(x) > 0){
        cres_gr <- cres_gr[-unique(x$queryHits)]
      } 
      
    }
  }
  
  if (exists("nbp") && !is.null(nbp)) {
    cat("Reading chromosome sizes...\n")
    chrom_sizes <- read.table(file = chrSizesFile, header = F, stringsAsFactors = F, sep = "\t")
    colnames(chrom_sizes) <- c("chr", "chr_size")
    
    # Adjust CREs
    cat("Adjusting CRE length...\n")
	
	# Need to work with genomic ranges object. The adjust_CREs function not suitable.
	# Instead have added GenomicRanges resize and trim functions that will complete this process.
	#cres <- adjust_CREs(cres, nbp, chrom_sizes)  
	
	cres_gr <- GenomicRanges::resize(cres_gr, width = nbp, fix = "center")

    
    idx <- which(chrom_sizes$chr %in% names(GenomeInfoDb::seqlengths(cres_gr)))
    if (length(idx) == length(names(GenomeInfoDb::seqlengths(cres_gr))))
    {
      suppressWarnings(GenomeInfoDb::seqlengths(cres_gr) <- chrom_sizes$chr_size[idx])
      cres_gr <- IRanges::trim(cres_gr)
	  # remove occurances that do not have a significant width
	  idx.small <- which(GenomicRanges::width(cres_gr) < nbp)
	  if (length(idx.small) > 0){
		  message(paste0("Removing ", length(idx.small), " regions as they were truncated to a length smaller than ", nbp, " nt (i.e. were located near the edge of a chromosome)"))
		  cres_gr <- cres_gr[idx.small * -1]
	  }
    } else {	
      warning("Cannot check out of bound granges as chromosome names between peaks and reference did not match")
    }
    
    
  }
  
  if(length(cres_gr) == 0){
    stop(paste("No remaining CREs after applying filters."))
  }	
  cres <- as.data.frame(cres_gr)[, c("seqnames","start", "end", "cellType")] 
  
  
  # Save filtered regions
  if (! is.null(out_bed)){
	  message(paste0("Saving ", length(cres_gr), " CREs...\n"))
	  write.table(x = cres, file = out_bed, quote = FALSE, col.names = FALSE
		      , row.names = FALSE, sep = "\t")
  }
  
  
  
  if (! is.null(celloutputDir)){
	  # Generate BED files for each cell type
	  # split into bed files
	  message("Preparing output directories and files so that everything is set for motif searching")
	  colnames(cres) <- c("chrom", "start", "end", "cellType")
    
    # Filter out low proportion of cells
    cellnumbers <- table(cres$cellType)
    idx.toRemove <- which(cellnumbers / sum(cellnumbers) <  (minCellPercent / 100) ) 
    idx.toRemove <- names(cellnumbers)[idx.toRemove]
	  
    if (length(idx.toRemove) > 0){
      cres <- cres[!cres$cellType %in% idx.toRemove, ]
      message(paste0("\nThe following cell(s) represent less than ", minCellPercent,
                     "% of all cells.", 
                     "\nThis is the predefined cutoff as defined by parameter minCellPercent.",
                     "\nOutputs for following cell(s) will therefore not be generated:\n", 
                     paste0(idx.toRemove, collapse = "\n"), "\n"))
    }
    
 
    
    for(i in unique(cres$cellType)){
	    fn <- paste0(celloutputDir, "/", gsub(" ", "_", i), ".bed")
	    write.table(x = cres[cres$cellType == i, ], file = fn, 
			quote = F, col.names = F, row.names = F, sep = "\t")
    }
    cellNames <- paste0(gsub(" ", "_", unique(cres$cellType)), ".bed", collapse = "\n")
    message(paste0("The following files have been prepared: \n", cellNames))
  }
  
  
}


#############################################################
#' textToBED
#' @description Converts a text file to BED format for peak data 
#'
#' @param inputTextFile  Input bed filename
#' @param header Boolean (default TRUE). Whether the input file has a header.
#' @param inputcolnames Column names corresponding to CRE coordinates and cell type/state annotation. Default to 'c("peak_chr","peak_start", "peak_end", "celltype_specificity")'.
#' @param sep Character separating columns. Default to ','.
#' @param outputFileName Output file name. Default to "out.bed".
#' @param removeMultiAnnotatedPeaks Boolean (default TRUE). Whether to remove CREs annotated to multiple cell types/states.
#' @param removeUnnanotatedPeaks Boolean (default TRUE). Whether to remove CREs not annotated to a cell type/state.
#' @examples 
#' \dontrun{
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' filename <- paste0(extdata_path,"/Pijuan_etal_table_S6.csv.gz")
#'
#' textToBed(inputTextFile = filename, outputFileName =  "mouseE8.25_peaks_filt.bed")
#' }
#'
#' @export
#'
textToBED <- function(inputTextFile = NULL, 
                      header = TRUE,
                      inputcolnames = c("peak_chr", "peak_start", "peak_end", "celltype_specificity"),
                      sep = ",",
                      outputFileName = "out.bed",
                      removeMultiAnnotatedPeaks = TRUE,
                      removeUnnanotatedPeaks = TRUE){
  cnameLength <- length(inputcolnames)
  if (cnameLength != 4){
	  error(paste0("Only ", cnameLength, " column names provided. Must be 4 column names defined.") )
  }
  
  txtData <- read.table(file = inputTextFile, header = header, sep = sep, stringsAsFactors = F)
  
  if (header == FALSE){
	  colnames(txtData)[1:4] <- inputcolnames
  }
  message("Input text file has ", nrow(txtData), " entries")
  
  
  
  if (removeUnnanotatedPeaks){
	  # removing all the peaks that were not annotated to a cell type
	  message("Removing ", length(which(is.na(txtData[,inputcolnames[4]]))), " entries that are not annotated with cell type")
	  txtData <- txtData[!is.na(txtData[, inputcolnames[4]]), ]  # celltype_specificity
  }
  
  if(removeMultiAnnotatedPeaks){
	  # remove peaks annotated to multiple cell types
	  idx <- paste(txtData[,inputcolnames[1]]
		       , paste(txtData[,inputcolnames[2]]
			       , txtData[,inputcolnames[3]], sep = "-")
		       , sep = ":")
    dup_idx <- (idx[duplicated(idx)])
    message("Removing ", length(dup_idx), " unique entries annotated to multiple cell types")
    remove_idx <- which(idx %in% dup_idx)
    txtData <- txtData[-remove_idx, ]
  }
  
  message("Creating bed file with ", nrow(txtData), " entries")
  
  # Save processed peak into a bed file
  write.table(x = txtData[, inputcolnames], file = outputFileName, col.names = F, row.names = F, quote = F, sep = "\t")
  
}

