# user_defined
#
# Retain reverse transcription stop/mutation/deletion information based on user-defined thresholds

#' Retain reverse transcription stop/mutation/deletion information based on user-defined thresholds
#' @param rocfile ROC input file of single sites information (with suffix _roc_plot.txt) generated from ground_truth() function
#' @param rRNAfile default we recommend hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
#' @param filtfile File to be filted
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
#'

user_defined <- function(rocfile, rRNAfile, filtfile, output_dir, output_name)
{
  #nohup Rscript $(dirname "$0")/user_defined_evaluation.r -f ${output_path}_roc_plot.txt -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed -t ${output_path}_user_defined_prediction_total.bed -o ${output_path} > ${output_path}_roc_bestthres.log 2>&1 &

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  outfile_prefix<-paste(output_dir,"/",output_name,sep="")

  # roc_plot.txt
  ROC_data<<-as.data.frame(fread(rocfile))
  colnames(ROC_data)<-c("chrom",
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
                        "extendSeq",
                        "pred_class",
                        "real_class")
  tab <- table(Predicted = ROC_data$pred_class,Actual = ROC_data$real_class)
  tab

  confusionMatrix(factor(ROC_data$pred_class,levels=c("0","1")),factor(ROC_data$real_class,levels=c("0","1")))
  ROC_data_sel<-ROC_data %>% select(treatPreRpmFold,treatAfterRpmFold,treatStopMeanFold,treatStopRatio,preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopMeanFoldRatio,real_class)
  print("total summary (rRNA psi and rRNA non-psi)")
  summary(ROC_data_sel)
  print("psi summary (all real rRNA psi)")
  real_rRNA_psi<-ROC_data_sel %>% filter(real_class=="1")
  summary(real_rRNA_psi)
  pdf(paste(outFile_prefix,"_real_rRNA_psi_datadensity.pdf",sep=""))
  datadensity(real_rRNA_psi, lwd = 1,group=cut2(real_rRNA_psi$treatStopRatio,g=2))#cut tretRtsRatio into 2 color group
  dev.off()

  #rRNA-psi-non-psi visualization
  ROC_data_melt<-melt(ROC_data[,c(11:38,41)],id.vars = "real_class")
  ROC_data_melt$real_class<-str_replace(as.character(ROC_data_melt$real_class),"0","non-psi")
  ROC_data_melt$real_class<-str_replace(as.character(ROC_data_melt$real_class),"1","psi")

  data_summary <- function(x) {
     m <- mean(x)
     ymin <- m-sd(x)
     ymax <- m+sd(x)
     return(c(y=m,ymin=ymin,ymax=ymax))
  }

  my_comparisons <- list( c("non-psi", "psi") )


  print("violin plotting...")
  tretRtsRatio_bp <- ggplot(ROC_data_melt%>%filter(str_detect(.$variable,"^tretRtsRatio$")), aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(tretRtsRatio)") +
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  rtsRatioFold_bp <- ggplot(ROC_data_melt%>%filter(str_detect(.$variable,"^rtsRatioFold$")), aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(rtsRatioFold)") +
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  tretBefRatio_bp <- ggplot(ROC_data_melt%>%filter(str_detect(.$variable,"^tretBefRatio$")), aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(tretBefRatio)")+
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  befRatioFold_bp <- ggplot(ROC_data_melt%>%filter(str_detect(.$variable,"^befRatioFold$")), aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(befRatioFold)")+
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  tretAftRatio_melt<-ROC_data_melt%>%filter(str_detect(.$variable,"^tretAftRatio$"))
  tretAftRatio_melt<-tretAftRatio_melt[tretAftRatio_melt$value!=0,]
  tretAftRatio_bp <- ggplot(tretAftRatio_melt, aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(tretAftRatio)") +
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  aftRatioFold_bp <- ggplot(ROC_data_melt%>%filter(str_detect(.$variable,"^aftRatioFold$")), aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(aftRatioFold)") +
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  tretMutRatio_bp <- ggplot(ROC_data_melt%>%filter(str_detect(.$variable,"^tretMutRatio$")), aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(tretMutRatio)")+
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  mutRatioFold_bp <- ggplot(ROC_data_melt%>%filter(str_detect(.$variable,"^mutRatioFold$")), aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(mutRatioFold)")+
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  tretDelRatio_bp <- ggplot(ROC_data_melt%>%filter(str_detect(.$variable,"^tretDelRatio$")), aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(tretDelRatio)") +
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  delRatioFold_bp <- ggplot(ROC_data_melt%>%filter(str_detect(.$variable,"^delRatioFold$")), aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(delRatioFold)")+
    scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))


  pdf(paste(outfile_prefix,"_ten_variables_rRNA_violinplot.pdf",sep=""),width=12,height=4)
  plot_grid(
    tretRtsRatio_bp,#'A'
    rtsRatioFold_bp,#'B'
    tretBefRatio_bp,#'C'
    befRatioFold_bp,#'D'
    tretAftRatio_bp,#'E'
    aftRatioFold_bp,#'F'
    tretMutRatio_bp,#'G'
    mutRatioFold_bp,#'H'
    tretDelRatio_bp,#'I'
    delRatioFold_bp,#'J'
    align = "hv",
    labels = c('A','B','C','D','E','F','G','H','I','J'),ncol=5,nrow=2)
  invisible(dev.off())

  #calculate evaluation indicators
  cat("\n\n=====================Calculate evaluation indicators=====================\n")
  confusion_matrix<-as.data.frame(tab)
  confusion_matrix$Predicted<-str_replace(confusion_matrix$Predicted,"1","psi")
  confusion_matrix$Predicted<-str_replace(confusion_matrix$Predicted,"0","non-psi")
  confusion_matrix$Actual<-str_replace(confusion_matrix$Actual,"1","psi")
  confusion_matrix$Actual<-str_replace(confusion_matrix$Actual,"0","non-psi")
  ud_TP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="psi",]$Freq#True Positives (TP)
  ud_FP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="non-psi",]$Freq#False Positives (FP)
  ud_TN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="non-psi",]$Freq#True Negatives (TN)
  ud_FN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="psi",]$Freq#False Negatives (FN)
  ud_TPR <- ud_TP / (ud_TP + ud_FN)#sensitivity (true positive rate, TPR)
  ud_TNR <- ud_TN / (ud_TN + ud_FP)#specifity (selectivity or true negative rate, TNR)
  ud_FPR <- 1-ud_TNR#False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
  ud_FNR <- 1-ud_TPR#False Negative Rate, FNR)
  ud_Prec <- ud_TP / (ud_TP + ud_FP)#Precision
  ud_Recall <- ud_TP / (ud_TP + ud_FN)#Recall
  ud_ACC <- (ud_TP + ud_TN) / (ud_TP + ud_TN + ud_FP + ud_FN)#accuracy
  ud_F1_score <- (2*ud_Recall*ud_Prec) / (ud_Recall + ud_Prec)#F1_score
  eval<-cbind(ud_TP,ud_FP,ud_TN,ud_FN,ud_TPR,ud_TNR,ud_FPR,ud_FNR,ud_Prec,ud_Recall,ud_ACC,ud_F1_score)
  eval<-round(eval,3)
  print(eval)
  write.table(eval,paste(outfile_prefix, '_ud_eval.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)

  #show ud evaluation as pdf table
  ud_eval_t_df<-as.data.frame(t(as.data.frame(eval)))
  colnames(ud_eval_t_df)<-"value_or_percentage"
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(fill = blues9[1:4], col=NA),
              fg_params=list(fontface=3)),
    colhead=list(fg_params=list(col="navyblue", fontface=4L)),
    rowhead=list(fg_params=list(col="orange", fontface=3L)))

  pdf(paste(outfile_prefix, '_user_defined_evaluation.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
  print(grid.arrange(
    tableGrob(ud_eval_t_df, theme=tt3),
    nrow=1))
  invisible(dev.off()) # Close the file

  #filt by best F1 score threshold
  to_filt<<-as.data.frame(fread(filtfile,skip=1))
  colnames(to_filt)<-c(
    "chrom",
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
    "extendSeq",
    "class")
  to_filt$pred_class<-str_replace(as.character(to_filt$pred_class),"0","non-psi")
  to_filt$pred_class<-str_replace(as.character(to_filt$pred_class),"1","psi")

  print(table(to_filt$pred_class))
  write.table(to_filt,paste(outfile_prefix, '_user_defined_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  final_pred<-to_filt %>% filter(pred_class=="psi")
  write.table(final_pred,paste(outfile_prefix, '_ud_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  write.table(final_pred,paste(outfile_prefix, '_ud_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

  #read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
  evidence<-read.table(rRNAfile,head=F)
  colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

  #final_pred miss/hit
  final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
  final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
  write.csv(final_pred_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiScan_ud_thres_hit.csv",sep=""))
  final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
  write.csv(final_pred_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiScan_ud_thres_miss.csv",sep=""))
  recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiScan+ud recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  print(paste("User-defined result in",output_dir,sep=" "))
}
