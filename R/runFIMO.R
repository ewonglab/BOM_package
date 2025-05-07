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
#' @param FIMO_update Binary variable to indicate if FIMO is updated. Default = FALSE.
#'
#' @examples
#'
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' bed_path <- paste0(extdata_path, "/tutorial/bed_files")
#'
#' # Replace /path/to/fimo with the path to FIMO tool
#' runFIMO(input_path = , motifs_path = paste0(extdata_path, "/gimme.vertebrate.v5.0.meme")
#' , out_path = "mouseE8.25_motifs", p_thresh = 0.0001
#' , FIMO_path = "/path/to/fimo")
#'
#' @export
#'
runFIMO <- function(input_path, motifs_path, out_path
                    , p_thresh = 0.0001, FIMO_path, FIMO_update = TRUE,  verbosity = 1){
  # create output directory
  mkidir_command <- paste("mkdir -p", out_path)
  system(mkidir_command)
  
  # get sequences for CRE coordinates
  message("Looking for motif instances with FIMO...")
  fastaFiles <- list.files(path = input_path, pattern = ".fa$", full.names = TRUE)
  # concise FIMO invocation
  pgc_flag <- if (as.logical(FIMO_update)) "--no-pgc" else ""
  fimo_exec <- file.path(FIMO_path, "fimo")

  for (fastaFile in fastaFiles) {
    out_dir <- file.path(out_path, sub("\\.fa$", "", basename(fastaFile)))
    args <- c(
      "--thresh",    p_thresh,
      "--verbosity", as.integer(verbosity),
      pgc_flag,
      "--o",         out_dir,
      motifs_path,
      fastaFile
    )
    cmd <- paste(
      shQuote(fimo_exec),
      paste(args[args != ""], collapse = " ")
    )
    system(cmd)
  }


                    }



