#############################################################
## runFIMO
#' Function to execute FIMO for at least two sets of CREs specific to different cell types/conditions 
#' 
#' @param input_path  Path to folder containing .bed files with CRE coordinates
#' @param genome_path  Path to genome file (fasta format)
#' @param motifs_path Motif database file (meme format)
#' @param out_path Output file folder
#' @param p_thresh p-value threshold (default 0.0001)
#' @param rm_fasta whether to remove the generated fasta files that contain CRE sequences (default TRUE) 
#'      
#' @examples 
#'
#' extdata_path <- system.file("extdata",package = "BagOfMotifs")
#' bed_path <- paste0(extdata_path, "/tutorial/bed_files")
#' 
#' 
#' runFIMO(input_path = bed_path, genome_path = "Mus_musculus_GRCm38.fa"
#' , motifs_path = "gimme.vertebrate.v5.0.meme", out_path = "./fimo_out")
#
#' @export

runFIMO <- function(input_path, genome_path, motifs_path, out_path
                    , p_thresh = 0.0001, rm_fasta = TRUE){
  # create output directory
  mkidir_command <- paste("mkdir -p", out_path)
  system(mkidir_command)
  
  # get sequences for CRE coordinates
  message("Extracting sequences for CRE coordinates and executing FIMO...\n")
  bed_files <- list.files(path = input_path, pattern = ".bed$", full.names = T)
  for(bed_file in bed_files){
    fa_file <- sub(".bed$", ".fa", bed_file)
    bedtools_command <- paste("bedtools getfasta -fi", genome_path
                              , "-bed", bed_file, "-fo", fa_file)
    system(bedtools_command)
    fimo_command <- paste("fimo --thresh", p_thresh, "--o"
                          , paste0(out_path, "/", sub(".bed$", ""
                                                      , basename(bed_file)))
                          , motifs_path, fa_file)
    system(fimo_command)
    if(rm_fasta){
      system(paste("rm", fa_file))
    }
  }
}
