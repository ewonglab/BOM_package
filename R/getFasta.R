#############################################################
## getFasta
#' Function to extract sequences for a set of CRE coordinates in BED format 
#' 
#' @param bed  Data frame containing CRE coordinates. The first three columns should correspond to chromosome, start and end positions. The fourth column should contain cell type/condition annotation.
#' @param genome  Genome as a BSgenome object
#' @param FastaFile Path to output fasta file. If NULL, a Biostrings object will be returned. Default = NULL. 
#'      
#' @examples 
#'
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' bed_path <- paste0(extdata_path, "/tutorial/bed_files/Cardiomyocytes.bed")
#' bed_cardiom <- read.table(file = bed_path, header = F, stringsAsFactors = F
#' , sep = '\t')  
#' 
#' # Match chromosome notation
#' bed$V1 <- paste0("chr", bed$V1)
#' 
#' # Load genome
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' Mmusculus <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' getFasta(bed = bed_cardiom, genome = Mmusculus, FastaFile = "Cardiomyocytes.fa")
#
#' @export
#' 
getFasta <- function(bed = NULL, genome, FastaFile = NULL){
  if(is.null(bed)){
    stop("CRE coordinates not provided.")
  }
  colnames(bed) <- c("chr", "start", "end")
  gRangesBed <-  with(bed, GRanges(chr, IRanges(start+1, end)))
  sequences <- Biostrings::getSeq(genome, gRangesBed)
  names(sequences) <- with(bed, paste(chr, paste(start, end, sep = "-"), sep = ":"))
  if(!is.null(FastaFile)){
    Biostrings::writeXStringSet(sequences, FastaFile)
  } else {
    return(sequences)
  }
}
