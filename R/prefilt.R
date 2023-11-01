# prefilt
#
# Filt out undesired psiScan information

#' Filt out undesired psiScan information
#' @param RpsiScan_res The R result of RpsiScan function
#' @param target_chrom The target chromsome to be retained
#' @param target_base The target base to be retained
#' @param tretRtsNum_thres The filtering threshold for tretRtsNum
#' @param tretBefNum_thres The filtering threshold for tretBefNum
#' @param tretAftNum_thres The filtering threshold for tretAftNum
#' @param tretMutNum_thres The filtering threshold for tretMutNum
#' @param tretDelNum_thres The filtering threshold for tretDelNum
#' @param output_dir The path to the output directory
#' @param output_name The output file name
#' @export
#'
prefilt <- function(RpsiScan_res,
                    target_chrom="^chr[0-9|a-z|A-Z]*$",
                    target_base="T",
                    tretRtsNum_thres=3,
                    tretBefNum_thres=3,
                    tretAftNum_thres=3,
                    tretMutNum_thres=3,
                    tretDelNum_thres=3,
                    output_dir,
                    output_name)
{
  prefilt_res<-RpsiScan_res %>% filter(grepl(target_chrom,`#chrom`),
                                       base==target_base,
                                       tretRtsNum > tretRtsNum_thres,
                                       tretBefNum > tretBefNum_thres,
                                       tretAftNum > tretAftNum_thres,
                                       tretMutNum > tretMutNum_thres,
                                       tretDelNum > tretDelNum_thres)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  output_dir_output_file<-paste(output_dir,"/",output_name,sep="")

  fwrite(x=prefilt_res,file=paste(output_dir_output_file,".txt",sep=""),sep="\t")
  write.table(prefilt_res,paste(output_dir_output_file,".bed",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  return(prefilt_res)
}
