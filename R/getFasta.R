#############################################################
## getFasta
#'
#' @param inputFile  input binary model predictions filename
#'
#' @examples 
#'
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' bed_path <- paste0(extdata_path, "/tutorial/bed_files")
#' bed <- list.files(path = bed_path, pattern = "(.*)bed$", full.names = F)[1]
#' Mmusculus <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' getFasta(bed, Mmusculus)
#
#' @export

getFasta <- function(bedFile, genome){
  bed <- read.table(file = bedFile, header = F, stringsAsFactors = F, sep = '\t')
  gRangesBed <-  with(bed, GRanges(V1, IRanges(V2+1, V3)))
  sequences <- Biostrings::getSeq(genome, gRangesBed)
  names(sequences) <- with(bed, paste(V1, paste(V2, V3, sep = "-"), sep = ":"))
  FastaFile <- paste0(sub("bed$", "", bedFile), "fa")
  Biostrings::writeXStringSet(sequences, FastaFile)
}
