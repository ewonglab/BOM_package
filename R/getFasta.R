#############################################################
## getFasta
#' Function to extract sequences for a set of CRE coordinates in BED format 
#'
#' @param bed  Data frame containing CRE coordinates. The first three columns should correspond to chromosome, start and end positions. 
#' The fourth column should contain cell type/condition annotation.
#' @param genome  Genome as a BSgenome object
#' @param FastaFile Path to output fasta file. If NULL, a Biostrings object will be returned. Default = NULL. 
#' @param UCSC  Logical. If TRUE, the chromosome names in the fasta file will be in UCSC format (e.g. chr1, chr2). If FALSE, the chromosome names will be in Ensembl format (e.g. 1, 2). Default = FALSE.
#'     
#' @examples 
#'
#' extdata_path <- system.file("extdata", package = "BagOfMotifs")
#' bed_path <- paste0(extdata_path, "/tutorial/bed_files/Cardiomyocytes.bed")
#' bed_cardiom <- read.table(file = bed_path, header = FALSE, stringsAsFactors = FALSE
#' , sep = "\t")  
#'
#' # Match chromosome notation
#' bed$V1 <- paste0("chr", bed$V1)
#'
#' # Load genome
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' Mmusculus <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' getFasta(bed = bed_cardiom, genome = Mmusculus, FastaFile = "Cardiomyocytes.fa")
#'
#' @export
#'
getFasta <- function(bed = NULL, genome, UCSC, FastaFile = NULL){
  if(is.null(bed)){
    stop("CRE coordinates not provided.")
  }
  colnames(bed) <- c("chr", "start", "end")
  gRangesBed <-  with(bed, GenomicRanges::GRanges(chr, IRanges::IRanges(start + 1, end)))
  GenomeInfoDb::seqlevelsStyle(gRangesBed) <- "ucsc"
  sequences <- Biostrings::getSeq(genome, gRangesBed)
  if (UCSC == TRUE){
	names(sequences) <- with(bed, paste(paste0("chr", chr), paste(start, end, sep = "-"), sep = ":"))
  }else{
	  names(sequences) <- with(bed, paste(chr, paste(start, end, sep = "-"), sep = ":"))
  }
  if(!is.null(FastaFile)){
    Biostrings::writeXStringSet(sequences, FastaFile)
  } else {
    return(sequences)
  }
}



#############################################################
## generateAllFasta
#' Wrapper for getFasta function. Will generate fasta files from bed files located within a specified directory
#' 
#' @param bedDir    Path to directory containing celltype/state specific BED files. The first three columns in the BED files should correspond to chromosome, start and end positions. The fourth column should contain cell type/state annotation.
#' @param genome    Genome as a BSgenome object
#' @param fastaDir  Directory to where fasta files are written to 
#' @param UCSC	 Logical. If TRUE, the chromosome names in the fasta file will be in UCSC format (e.g. chr1, chr2). If FALSE, the chromosome names will be in Ensembl format (e.g. 1, 2). Default = FALSE.
#'
#' @examples 
#'
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' bed_path <- paste0(extdata_path, "/tutorial/bed_files/")
#' fasta_path <- paste0(extdata_path, "/tutorial/fasta_files/")
#'
#' # Load genome
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' Mmusculus <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'
#' generateAllFasta(bedDir = bed_path, genome = Mmusculus, fastaDir = fasta_path)
#'
#' @export
#'
generateAllFasta <- function(bedDir = NULL, genome = NULL, fastaDir = NULL, UCSC = FALSE){
	if(is.null(bedDir)){
		stop("Bed input directory not provided.")
	}
	if(is.null(genome)){
		stop("Genome object not provided.")
	}
	if(is.null(fastaDir)){
		stop("Fasta output directory not provided.")
	}

	fl <- list.files(path = bedDir, pattern = "*.bed$")

	if (length(fl) == 0){
		stop(paste0("No Bed files found in following directory:\n",bedDir))
	}
	
	if (! dir.exists(fastaDir)){
		message("Attempting to create fasta output directory")
		dir.create(fastaDir)
	}
	
	allfastaFiles <- gsub(pattern = "\\.bed$", replacement = ".fa", x = fl)
	message(paste0("Creating ", length(allfastaFiles)," output fasta files"))
	for(i in 1:length(allfastaFiles)){
		bed_data <- read.table(file = paste0(bedDir, "/", fl[i]), header = F, stringsAsFactors = F, sep = "\t")  
		getFasta(bed = bed_data, genome, UCSC, FastaFile = paste0(fastaDir, "/", allfastaFiles[i]))
	}
	message("Fasta files generated")
}
