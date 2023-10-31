# prefilt
#
# Filt out undesired psiScan information

#' Filt out undesired psiScan information
#' @param RpsiScan_res The R result of RpsiScan function
#' @param target_chrom The target chromsome to be retained
#' @param target_base The target base to be retained
#' @param output_dir The path to the output directory
#' @param output_name The output file name
#' @export
#'
prefilt <- function(RpsiScan_res, target_chrom="^chr[0-9|a-z|A-Z]*$", target_base="T",output_dir,output_name)
{
  #awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $7=="T"){print $0}}' $out > ${outFileName}_roc_filt.bed
  prefilt_res<-RpsiScan_res %>% filter(grepl(target_chrom,`#chrom`),base==target_base)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  output_dir_output_file<-paste(output_dir,"/",output_name,sep="")

  fwrite(x=prefilt_res,file=paste(output_dir_output_file,".txt",sep=""),sep="\t")
  write.table(prefilt_res,paste(output_dir_output_file,".bed",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  return(prefilt_res)
}
