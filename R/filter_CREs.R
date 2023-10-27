############################################################
### adjust_CREs
#' Function to trim (or extend) CREs 
#' @param x
#' @param N
#'
#' @return  
adjust_CREs <- function(x, N){
  
  x$width <- with(x, V3-V2)
  #x$centre <- with(x, V2 + (width/2))
  x$centre <- with(x, V2 + round(width/2))
  #x$start <- ceiling(x$centre - N/2)
  x$start <- round(x$centre - N/2)
  #x$end <- ceiling(x$centre + N/2)
  x$end <- round(x$centre + N/2)

  # add chromosome size
  x$ord <- 1:nrow(x)
  x <- merge(x, chrom_sizes, by.x="V1", by.y="chr")

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
  
  #return to original order
  x <- x[with(x, order(ord)), ]
  x <- x[,colnames(x)[!colnames(x) %in% c("chr_size", "centre", "width", "V2", "V3", "ord")]]
  x <- x[,c("V1", "start", "end", setdiff(colnames(x), c("V1", "start", "end")))]
  return(x)
  
}




#############################################################
## filterCREs
#'
#' @param inputBedFile  input bed filename
#' @param annotFile     annotation file name
#' @param u number of basepairs upstream of transcriptional start sites
#' @param d number of basepairs downstream of transcriptional start site
#' @param nbp  number of base pairs
#' @param keep_proximal  boolean. If True keeps proximal regions close to Transcription start site
#' @param remove_proximal 
#' @param non_exonic
#' @param out_bed   name of output bed file
#' @returns
#'
#' @examples 
#'
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' input.bed <- paste0(extdata_path, "/tutorial/mouseE8.25_peaks.bed")
#' annot.file <- paste0(extdata_path,"/Mus_musculus.GRCm38.92.gtf.gz")
#' mouseChrSizes <- paste0(extdata_path,"/Mus_musculus.GRCm38.92.gtf.gz")
#' 
#' 
#'  filterCRE(inputBedFile = input.bed, annotFile = annot.file, 
#'                     remove_proximal = TRUE,  non_exonic = TRUE,
#' 					   u =  1000, d=1000, nbp=500,
#'                     chrSizes = mouseChrSizes, 
#'                     out_bed =  'mouseE8.25_peaks_filt.bed')
#
#' @export
filterCREs <- function(input_bed = NULL, 
				annot = NULL,
				chr_sizes = NULL,
				u = NULL,
				d = NULL,
				nbp = NULL,
				keep_proximal = FALSE,
				remove_proximal = FALSE,
				non_exonic = FALSE,
				out_bed = NULL)
{
	message("Reading CREs...\n")
	cres <- read.table(file = input_bed, header = F, stringsAsFactors = F, sep = '\t')
	cres_gr <- with(cres, GenomicRanges::GRanges(V1, IRanges::IRanges(V2+1, V3)))

	if(keep_proximal | remove_proximal | non_exonic){

	  message("Reading genome annotation...\n")
	  txdb_obj <- makeTxDbFromGFF(file = annot, format = "gtf")

	  if(non_exonic){
		message("Removing exonic regions...\n")
		exons <- exons(txdb_obj)
		x <- as.data.frame(findOverlaps(cres_gr, exons))
		if(nrow(x) > 0){
		  cres <- cres[-unique(x$queryHits),]
		} 
	  }
	  
	  if(keep_proximal & remove_proximal){
		stop(paste("'keep_proximal' and 'remove_proximal' are mutually exclusive"))
	  }
	  if(keep_proximal){
		message("Keeping proximal regions to TSSs...\n")
		proximal <- promoters(x = txdb_obj, upstream = u, downstream = d)
		cres_gr <- with(cres, GRanges(V1, IRanges(V2+1, V3)))
		x <- as.data.frame(findOverlaps(cres_gr, proximal))
		if(nrow(x) > 0){
		  cres <- cres[unique(x$queryHits),]
		} else {
		  stop(paste("No remaining CREs"))
		}
	  }
	  if(remove_proximal){
		message("Removing proximal regions to TSSs...\n")
		proximal <- promoters(x = txdb_obj, upstream = u, downstream = d)
		cres_gr <- with(cres, GRanges(V1, IRanges(V2+1, V3)))
		x <- as.data.frame(findOverlaps(cres_gr, proximal))
		if(nrow(x) > 0){
		  cres <- cres[-unique(x$queryHits),]
		} 
		
	  }
	}

	#if((u == NULL & d != NULL) | (d == NULL & u != NULL)){
	#  stop(paste("u AND d should be provided to adjust CREs"))
	if (exists("nbp") && !is.null(nbp)) {
	  cat("Reading chromosome sizes...\n")
	  chrom_sizes <- read.table(file = chr_sizes, header = F, stringsAsFactors = F, sep ='\t')
	  colnames(chrom_sizes) <- c("chr", "chr_size")
	  
	  # Adjust CREs
	  cat("Adjusting CRE length...\n")
	  cres <- adjust_CREs(cres, nbp)
	}

	# Save filtered regions
	message("Saving CREs...\n")

	if(nrow(cres) == 0){
	  stop(paste("No remaining CREs"))
	}

	write.table(x = cres, file = out_bed, quote = F, col.names = F
				, row.names = F, sep ='\t')


}






