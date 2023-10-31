# RpsiACAscan
#
# This function run psiACAscan and return systematic H/ACA snoRNA-target RNA interaction information

#' This function run psiACAscan and return systematic H/ACA snoRNA-target RNA interaction information
#' @param acaboxseq The H/ACA box sequence file
#' @param modifiedfile The Modified file to be scanned (provided by mod_info() function)
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
RpsiACAscan <- function(acaboxseq,
                        modifiedfile,
                        output_dir,
                        output_name)
{
  psiACAscan_path <- system.file("program", "psiACAscan", package = "RpsiScan")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  output_dir_output_name<-paste(output_dir,"/",output_name,sep="")

  paste(output_dir_output_name,"_modinfo.txt",sep="")

  cmd <- paste(psiACAscan_path,
               "-f", acaboxseq,
               "-m", modifiedfile,
               ">", output_dir_output_name,
               sep=" ")

  print(cmd)
  system(cmd)
  # assign(output_name,fread(output_dir_output_name),envir = .GlobalEnv)
}
