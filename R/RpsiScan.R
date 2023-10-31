# RpsiScan
#
# This function run psiScan and return overall reverse transcription stop/mutation/deletion information

#' This function run psiScan and return overall reverse transcription stop/mutation/deletion information
#' @param treat_bam The path to the treat bam file
#' @param input_bam The path to the input bam file
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @param options Additional options to pass to the program
#' @export
#'
RpsiScan <- function(genome_fa,
                       genome_fai,
                       treat_bam,
                       input_bam,
                       output_dir,
                       output_name,
                       options = "-p 1.5 -t 5 -r 0.05 -M 1 -f 1 -s -n -w 20")
{
  psiScan_path <- system.file("program", "psiScan", package = "RpsiScan")
  system(paste("chmod 777", psiScan_path, sep=" "))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  output_dir_output_name<-paste(output_dir,"/",output_name,sep="")
  cmd <- paste(psiScan_path,
               "--fa", genome_fa,
               "--fai", genome_fai,
               "--treat", treat_bam,
               "--input", input_bam,
               options,
               "-o", paste(output_dir_output_name,".txt",sep=""),
               paste("2>", output_dir, "/", output_name,".log",sep=""),
               sep=" ")

  print(cmd)
  system(cmd)
  assign(output_name,fread(paste(output_dir_output_name,".txt",sep="")),envir = .GlobalEnv)
}
