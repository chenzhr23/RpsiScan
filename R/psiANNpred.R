# psiANNpred
#
# Build ANN model based on ground truth pseudouridylation data set

#' Build ANN model based on ground truth pseudouridylation data set
#' @param filtfile File to be filted (File ready to be predicted)
#' @param annmodelfile ANN model file
#' @param rRNAfile default we recommend hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
#'

psiANNpred <- function(filtfile, annmodelfile, rRNAfile, output_dir,output_name)
{
  #nohup Rscript $(dirname "$0")/ann_predict_totalRNA.r -f ${input%.bed}_ann_filt.bed -k $ANN_model -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -o ${input%.bed} > ${input%.bed}.log

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  outfile_prefix<-paste(output_dir,"/",output_name,sep="")

  load(file = annmodelfile)#"my-ann.RData"

  #filt by ann best model
  cat("\n\n=====================Filt by ann best model=====================\n")
  cat("below is your input data ready to be predicted...\n")
  # get prediction
  to_pred<-as.data.frame(fread(filtfile))
  colnames(to_pred)<-c("chrom",
                       "chromStart",
                       "chromEnd",
                       "name",
                       "foldChange",
                       "strand",
                       "base",
                       "rtsPval",
                       "mutPval",
                       "delPval",
                       "tretRtsNum",
                       "tretRtsRpm",
                       "ctrlRtsNum",
                       "ctrlRtsRpm",
                       "tretHgtNum",
                       "ctrlHgtNum",
                       "tretBefNum",
                       "ctrlBefNum",
                       "tretAftNum",
                       "ctrlAftNum",
                       "rpmFold",
                       "tretRtsRatio", #*
                       "ctrlRtsRatio",
                       "rtsRatioFold", #*
                       "tretBefRatio", #*
                       "ctrlBefRatio",
                       "befRatioFold", #*
                       "tretAftRatio", #*
                       "ctrlAftRatio",
                       "aftRatioFold", #*
                       "tretMutNum",
                       "ctrlMutNum",
                       "tretMutRatio", #*
                       "ctrlMutRatio",
                       "mutRatioFold", #*
                       "delPos",
                       "tretDelNum",
                       "ctrlDelNum",
                       "tretDelRatio", #*
                       "ctrlDelRatio",
                       "delRatioFold", #*
                       "extendSeq")
  to_pred$base<-str_replace_all(to_pred$base, "TRUE", "T")
  to_pred_var<-to_pred %>% select(tretRtsRatio,rtsRatioFold,tretBefRatio,befRatioFold,tretAftRatio,aftRatioFold,tretMutRatio,mutRatioFold,tretDelRatio,delRatioFold)
  to_pred_var = scale(to_pred_var, center = training_set_mean, scale = training_set_sd)
  str(to_pred_var)

  # Predicting Result
  ann_data_sel.prediction <- neuralnet::compute(ANN_data_sel.net, to_pred_var)
  idx <- apply(ann_data_sel.prediction$net.result, 1, which.max)
  net.result <- as.data.frame(ann_data_sel.prediction$net.result)
  colnames(net.result) <- c('psi', 'non_psi')
  pred <- c('psi', 'non_psi')[idx]
  table(pred)

  #add attr to ann_pred_data
  ann_pred_data <- to_pred
  ann_pred_data <- cbind(ann_pred_data, net.result, pred)
  write.xlsx(ann_pred_data, paste(outfile_prefix, '_ann_total_prediction.xlsx', sep = ""), overwrite = TRUE)
  write.table(ann_pred_data,paste(outfile_prefix, '_ann_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  write.table(ann_pred_data,paste(outfile_prefix, '_ann_total_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

  final_pred <- ann_pred_data[ann_pred_data$pred == "psi",]
  print("psi predicted by nn model and fold change (foldChange) >2:")
  final_pred <- final_pred[final_pred$foldChange > 2,]
  print("Final psi prediction:")
  table(final_pred$pred)
  write.xlsx(final_pred, paste(outfile_prefix, '_ann_psi_prediction.xlsx', sep = ""), overwrite = TRUE)
  write.table(final_pred,paste(outfile_prefix, '_ann_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  write.table(final_pred,paste(outfile_prefix, '_ann_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

  #known evidence
  evidence <- read.table(rRNAfile, head = F)
  colnames(evidence) <- c("chrom", "chromStart", "chromEnd", "rRNA_anno", "score", "strand")
  evidence$rRNA_uniq_id <- paste(evidence$chrom, evidence$chromStart, evidence$chromEnd, evidence$strand, sep = "_")
  final_pred$rRNA_uniq_id <- paste(final_pred$chrom, final_pred$chromStart, final_pred$chromEnd, final_pred$strand, sep = "_")
  final_pred_evidence <- final_pred %>% left_join(evidence, by = c("rRNA_uniq_id" = "rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
  write.csv(final_pred_evidence, paste(outfile_prefix, "_rtsSeeker_ann_hit.csv", sep = ""))
  final_pred_no_evidence <- evidence %>% left_join(final_pred, by = c("rRNA_uniq_id" = "rRNA_uniq_id")) %>% filter(is.na(base))
  write.csv(final_pred_no_evidence, paste(outfile_prefix, "_rtsSeeker_ann_miss.csv", sep = ""))
  recovery <- paste(round(length(unique(final_pred_evidence$rRNA_uniq_id)) / length(unique(evidence$rRNA_uniq_id)) * 100, 2), "%", sep = "")
  cat("rtsSeeker+ANN recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)", recovery, "rRNA psi sites in all known chrom21\n")

  print(paste("ANN prediction in",output_dir,sep=" "))
}
