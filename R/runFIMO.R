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
