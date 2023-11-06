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


#############################################################
## runFIMO
#' Function to execute FIMO for at least two sets of CREs specific to different cell types/conditions 
#' 
#' @param input_path  Path to folder containing fasta (.fa) files with CRE sequences
#' @param motifs_path Motif database file (meme format)
#' @param out_path Output file folder. The output folder will be created and FIMO output for each fasta file will be saved to this directory.
#' @param p_thresh p-value threshold (default 0.0001)
#' @param FIMO_path Path to FIMO
#' @param verbosity Level of verbosity for FIMO. Default = 1.
#'      
#' @examples 
#'
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' bed_path <- paste0(extdata_path, "/tutorial/bed_files")
#' 
#' # Replace /path/to/fimo with the path to FIMO tool
#' runFIMO(input_path = , motifs_path = "gimme.vertebrate.v5.0.meme"
#' , out_path = "mouseE8.25_motifs", p_thresh = 0.0001
#' , FIMO_path = '/path/to/fimo')
#
#' @export
#' 
runFIMO <- function(input_path, motifs_path, out_path
                    , p_thresh = 0.0001, FIMO_path, verbosity = 1){
  # create output directory
  mkidir_command <- paste("mkdir -p", out_path)
  system(mkidir_command)
  
  # get sequences for CRE coordinates
  message("Looking for motif instances with FIMO...")
  fastaFiles <- list.files(path = input_path, pattern = ".fa$", full.names = T)
  for(fastaFile in fastaFiles){
    fimo_command <- paste(paste(FIMO_path, "fimo", sep = '/')
                          , " --thresh", p_thresh, "--verbosity", as.integer(verbosity)
                          , "--o", paste0(out_path, "/",
                                          sub(".fa$", "", basename(fastaFile)))
                          , motifs_path, fastaFile)
    system(fimo_command)
  }
}
