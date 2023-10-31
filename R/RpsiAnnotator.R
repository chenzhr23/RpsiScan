# RpsiAnnotator
#
# This function run psiAnnotator and add annotation for reverse transcription stop information

#' This function run psiAnnotator and add annotation for reverse transcription stop information
#' @param input_bedfile The input bed file to be annotated
#' @param annotation_bed6file The bed6 file for annotation
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
#'

add_biotype<-function(anno_bed_file,output_name){
  search_list<-list(
    data.frame(Gene_biotype="mRNA",Gene_biotype_detail=c("protein_coding")),
    data.frame(Gene_biotype="pseudogene",Gene_biotype_detail=c("rRNA_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_processed_pseudogene", "unitary_pseudogene", "pseudogene", "polymorphic_pseudogene", "transcribed_unitary_pseudogene", "TR_V_pseudogene", "TR_J_pseudogene", "IG_V_pseudogene", "IG_C_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "IG_J_pseudogene")),
    data.frame(Gene_biotype="lncRNA",Gene_biotype_detail=c("lncRNA","processed_transcript","lincRNA","non_coding","3prime_overlapping_ncRNA","3prime_overlapping_ncrna","sense_intronic","antisense","sense_overlapping","known_ncrna","macro_lncRNA","bidirectional_promoter_lncRNA","retained_intron","TEC")),
    data.frame(Gene_biotype="sncRNA",Gene_biotype_detail=c("snRNA","snoRNA","misc_RNA","miscRNA","miRNA","ribozyme","sRNA","scRNA","scaRNA","srpRNA","tRNA-Deu","tRNA-RTE","piRNA","siRNA")),
    data.frame(Gene_biotype="rRNA",Gene_biotype_detail=c("rRNA","Mt_rRNA")),
    data.frame(Gene_biotype="tRNA",Gene_biotype_detail=c("tRNA","Mt_tRNA","vaultRNA")),
    data.frame(Gene_biotype="IG_gene",Gene_biotype_detail=c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene")),
    data.frame(Gene_biotype="TR_gene",Gene_biotype_detail=c("TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene")),
    data.frame(Gene_biotype="repeatMasker",Gene_biotype_detail=c("5S-Deu-L2","Alu","centr","CR1","DNA","DNA?","ERV1","ERV1?","ERVK","ERVL","ERVL?","ERVL-MaLR","Gypsy","Gypsy?","hAT","hAT?","hAT-Ac","hAT-Blackjack","hAT-Charlie","hAT-Tag1","hAT-Tip100","hAT-Tip100?","Helitron","Helitron?","L1","L2","Low_complexity","LTR","LTR?","MIR","MULE-MuDR","nonsense_mediated_decay","non_stop_decay","Penelope","PiggyBac","PiggyBac?","RNA","RTE-BovB","RTE-X","Satellite","Simple_repeat","SVA","TcMar?","TcMar-Mariner","TcMar-Tc2","TcMar-Tigger","telo","Unknown","acro","Crypton","Dong-R4","I-Jockey","Kolobok","L1-Tx1","Merlin","MULE-MuDR?","PIF-Harbinger","SINE?","TcMar","TcMar-Pogo","TcMar-Tc1")),
    data.frame(Gene_biotype="intergenic",Gene_biotype_detail=c("intergenic")),
    data.frame(Gene_biotype="circRNA",Gene_biotype_detail=c("circRNA"))
  )

  search_df<-do.call(rbind,search_list)

  anno_bed_file_input<-as.data.frame(fread(anno_bed_file,fill=TRUE))
  anno_bed_file_input_separate<-anno_bed_file_input %>%
    separate_wider_delim(V7, "|", names = c("Transcript_ID",
                                            "Transcript_NAME",
                                            "Gene_ID",
                                            "Gene_NAME",
                                            "Gene_Type",
                                            "Region"),too_few = "align_start")

  anno.biotype.bed<-anno_bed_file_input_separate %>% left_join(search_df,by=c("Gene_Type"="Gene_biotype_detail"))
  anno.biotype.bed$V8<-str_replace_all(anno.biotype.bed$V8,"0","intergenic")

  anno_bed_file_Gene_biotype<-sub("_anno.bed","_anno.biotype.bed",anno_bed_file)
  fwrite(x=anno.biotype.bed,file = anno_bed_file_Gene_biotype,sep="\t",col.names = FALSE)
  assign(paste(output_name,"_anno_biotype_bed_R",sep=""),fread(anno_bed_file_Gene_biotype),envir = .GlobalEnv)
  return(anno_bed_file_Gene_biotype)

}

remove_redundancy<-function(anno.biotype.bed.file, input.bed.file, outfile){

  # anno.biotype.bed<-read.table(anno.biotype.bed.file,sep="\t")
  anno.biotype.bed<-as.data.frame(fread(anno.biotype.bed.file))
  anno.biotype.bed_uniqid<-paste(anno.biotype.bed$V1,anno.biotype.bed$V2,anno.biotype.bed$V3,anno.biotype.bed$V6,sep="_")
  # input.bed<-read.table(input.bed.file,sep="\t")
  input.bed<-as.data.frame(fread(input.bed.file))
  input.bed_uniqid<-paste(input.bed$V1,input.bed$V2,input.bed$V3,input.bed$V6,sep="_")
  probe<-which(input.bed_uniqid %in% anno.biotype.bed_uniqid)
  input.bed<-input.bed[probe,]
  input.bed$input.bed_uniqid<-input.bed_uniqid
  anno.biotype.bed$anno.biotype.bed_uniqid<-anno.biotype.bed_uniqid
  add_seq <- anno.biotype.bed %>% left_join(input.bed,by=c("anno.biotype.bed_uniqid"="input.bed_uniqid"))
  add_seq<-add_seq %>% select(-anno.biotype.bed_uniqid)
  seq_index <- which(as.data.frame(unlist(apply(add_seq,2,function(x){unique(str_detect(x[1],"^[AGCT].*[AGCT]$")&str_length(x)>1)})))=="TRUE")


  extendSeq_split<-as.data.frame(str_split_fixed(add_seq[,seq_index], "Y", 2))
  colnames(extendSeq_split)<-c("extendSeq_bef","extendSeq_aft")
  add_seq<-cbind(add_seq,extendSeq_split)
  seq_group<-data.frame(uniq_seq=unique(add_seq$extendSeq_aft))
  seq_group$seq_id<-1:length(seq_group$uniq_seq)
  add_seq <- add_seq %>% left_join(seq_group,by=c("extendSeq_aft"="uniq_seq"))


  priority<-c("rRNA","tRNA","sncRNA","mRNA","lncRNA","circRNA","IG_gene","TR_gene","pseudogene","repeatMasker","intergenic")
  priority_table<-data.frame(type=priority,priority_rank=c(seq(length(priority))))
  add_seq <- add_seq %>% left_join(priority_table,by=c("V19.x"="type"))
  add_seq$group_id<-paste("group",add_seq$seq_id,add_seq$priority_rank,sep="-")
  write.table(add_seq,paste(outfile,"_add_seq_group.bed",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  # write.xlsx(add_seq,paste(outfile,"_add_seq_group.xlsx",sep=""), overwrite = TRUE)

  add_seq_test<<-add_seq
  add_seq_group_list<-split(add_seq,add_seq$seq_id)

  chro_id<-c("chr21",paste("chr",1:20,sep=""),"chr22","chrM","chrX","chrY")
  chro_priority<-data.frame(chro_id=chro_id,chro_rank=1:length(chro_id))
  invisible(lapply(seq_along(add_seq_group_list),function(x){
    add_seq_group_list_tmp <- add_seq_group_list[[x]] %>% left_join( chro_priority,by=c("V1.x"="chro_id"))
    add_seq_group_list_tmp<-arrange(add_seq_group_list_tmp,priority_rank,chro_rank)
    add_seq_group_list_tmp<-add_seq_group_list_tmp %>% select(-chro_rank)
    add_seq_group_list[[x]]<<-head(add_seq_group_list_tmp,1)
  }))

  add_seq_group_uniq<<-do.call(rbind.data.frame, add_seq_group_list)

  if(str_detect(input.bed.file,"roc_psi_prediction.bed")){
    colnames(add_seq_group_uniq)<-c("chrom",
                                    "chromStart",
                                    "chromEnd",
                                    "name",
                                    "foldChange",
                                    "strand",
                                    "Transcript_ID",
                                    "Transcript_NAME",
                                    "Gene_ID",
                                    "Gene_NAME",
                                    "Gene_Biotype_detail",
                                    "Region",
                                    "Gene_Feature",
                                    "tolBaseNum",
                                    "tqueryCov",
                                    "tsampCov",
                                    "tupDist",
                                    "tdownDist",
                                    "Gene_Biotype",
                                    "chrom.y",
                                    "chromStart.y",
                                    "chromEnd.y",
                                    "name.y",
                                    "foldChange.y",
                                    "strand.y",
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
                                    "roc_pred",
                                    "extendSeq_bef",
                                    "extendSeq_aft",
                                    "seq_id",
                                    "priority_rank",
                                    "group_id")
    add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")
  }else if(str_detect(input.bed.file,"svm_psi_prediction.bed")){
    colnames(add_seq_group_uniq)<-c("chrom",
                                    "chromStart",
                                    "chromEnd",
                                    "name",
                                    "foldChange",
                                    "strand",
                                    "Transcript_ID",
                                    "Transcript_NAME",
                                    "Gene_ID",
                                    "Gene_NAME",
                                    "Gene_Biotype_detail",
                                    "Region",
                                    "Gene_Feature",
                                    "tolBaseNum",
                                    "tqueryCov",
                                    "tsampCov",
                                    "tupDist",
                                    "tdownDist",
                                    "Gene_Biotype",
                                    "chrom.y",
                                    "chromStart.y",
                                    "chromEnd.y",
                                    "name.y",
                                    "foldChange.y",
                                    "strand.y",
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
                                    "svm_pred.decision.values",
                                    "svm_pred_desc_class",
                                    "svm_psi_prob",
                                    "svm_non_psi_prob",
                                    "svm_pred",
                                    "svm_pred_prob_class",
                                    "extendSeq_bef",
                                    "extendSeq_aft",
                                    "seq_id",
                                    "priority_rank",
                                    "group_id")
    add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")
  }else if(str_detect(input.bed.file,"ann_psi_prediction.bed")){
    colnames(add_seq_group_uniq)<-c("chrom",
                                    "chromStart",
                                    "chromEnd",
                                    "name",
                                    "foldChange",
                                    "strand",
                                    "Transcript_ID",
                                    "Transcript_NAME",
                                    "Gene_ID",
                                    "Gene_NAME",
                                    "Gene_Biotype_detail",
                                    "Region",
                                    "Gene_Feature",
                                    "tolBaseNum",
                                    "tqueryCov",
                                    "tsampCov",
                                    "tupDist",
                                    "tdownDist",
                                    "Gene_Biotype",
                                    "chrom.y",
                                    "chromStart.y",
                                    "chromEnd.y",
                                    "name.y",
                                    "foldChange.y",
                                    "strand.y",
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
                                    "ann_psi_prob",
                                    "ann_non_psi_prob",
                                    "ann_pred",
                                    "extendSeq_bef",
                                    "extendSeq_aft",
                                    "seq_id",
                                    "priority_rank",
                                    "group_id")
    add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")
  }else if (str_detect(input.bed.file,"user_defined_psi_prediction.bed")){
    colnames(add_seq_group_uniq)<-c("chrom",
                                    "chromStart",
                                    "chromEnd",
                                    "name",
                                    "foldChange",
                                    "strand",
                                    "Transcript_ID",
                                    "Transcript_NAME",
                                    "Gene_ID",
                                    "Gene_NAME",
                                    "Gene_Biotype_detail",
                                    "Region",
                                    "Gene_Feature",
                                    "tolBaseNum",
                                    "tqueryCov",
                                    "tsampCov",
                                    "tupDist",
                                    "tdownDist",
                                    "Gene_Biotype",
                                    "chrom.y",
                                    "chromStart.y",
                                    "chromEnd.y",
                                    "name.y",
                                    "foldChange.y",
                                    "strand.y",
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
                                    "ud_pred",
                                    "extendSeq_bef",
                                    "extendSeq_aft",
                                    "seq_id",
                                    "priority_rank",
                                    "group_id")
    add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")
  }else{
    colnames(add_seq_group_uniq)<-c("chrom",
                                    "chromStart",
                                    "chromEnd",
                                    "name",
                                    "foldChange",
                                    "strand",
                                    "Transcript_ID",
                                    "Transcript_NAME",
                                    "Gene_ID",
                                    "Gene_NAME",
                                    "Gene_Biotype_detail",
                                    "Region",
                                    "Gene_Feature",
                                    "tolBaseNum",
                                    "tqueryCov",
                                    "tsampCov",
                                    "tupDist",
                                    "tdownDist",
                                    "Gene_Biotype",
                                    "chrom.y",
                                    "chromStart.y",
                                    "chromEnd.y",
                                    "name.y",
                                    "foldChange.y",
                                    "strand.y",
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
                                    "extendSeq_bef",
                                    "extendSeq_aft",
                                    "seq_id",
                                    "priority_rank",
                                    "group_id")
    add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")
  }

  if(dim(add_seq_group_uniq)[1]>0){
    pseudoU_anno_genetype_num<-arrange(as.data.frame(table(add_seq_group_uniq$Gene_Biotype)),desc(Freq))
    write.table(pseudoU_anno_genetype_num,paste(outfile,"_anno_gene_biotype_num.sort",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
    percent<-round(100*pseudoU_anno_genetype_num$Freq/sum(pseudoU_anno_genetype_num$Freq),2)
    percent <-paste('(',percent, "%", ", ", pseudoU_anno_genetype_num$Freq,')', sep = "")
    pseudoU_anno_genetype_num$Var1 <- paste(pseudoU_anno_genetype_num$Var1, percent, sep = '')
    pseudoU_anno_genetype_num$Var1 <- factor(pseudoU_anno_genetype_num$Var1, levels = pseudoU_anno_genetype_num$Var1)
    mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"),brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Pastel1"),brewer.pal(8, "Pastel2"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))
    pie_plot<-ggplot(data = pseudoU_anno_genetype_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) +
      geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +
      theme(axis.text = element_blank()) +
      theme(axis.ticks = element_blank()) +
      theme(panel.grid = element_blank(),
            panel.background=element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank()) +
      guides(fill = guide_legend(title = "Gene Biotype"))+
      scale_fill_manual(values = mycol)

    pdf(paste(outfile,"_gene_biotype_piechart.pdf",sep=""))
    print(pie_plot)
    invisible(dev.off())

    pseudoU_anno_genetype_num<-arrange(as.data.frame(table(add_seq_group_uniq$Gene_Feature)),desc(Freq))
    percent<-round(100*pseudoU_anno_genetype_num$Freq/sum(pseudoU_anno_genetype_num$Freq),2)
    percent <-paste('(',percent, "%", ", ", pseudoU_anno_genetype_num$Freq,')', sep = "")
    pseudoU_anno_genetype_num$Var1 <- paste(pseudoU_anno_genetype_num$Var1, percent, sep = '')
    pseudoU_anno_genetype_num$Var1 <- factor(pseudoU_anno_genetype_num$Var1, levels = pseudoU_anno_genetype_num$Var1)
    mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"),brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Pastel1"),brewer.pal(8, "Pastel2"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))
    pie_plot<-ggplot(data = pseudoU_anno_genetype_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) +
      geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +
      theme(axis.text = element_blank()) +
      theme(axis.ticks = element_blank()) +
      theme(panel.grid = element_blank(),
            panel.background=element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank()) +
      guides(fill = guide_legend(title = "Gene Feature"))+
      scale_fill_manual(values = mycol)

    pdf(paste(outfile,"_gene_feature_piechart.pdf",sep=""))
    print(pie_plot)
    invisible(dev.off())


    write.table(add_seq_group_uniq,paste(outfile,"_add_seq_group_uniq.bed",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
    write.xlsx(add_seq_group_uniq,paste(outfile,"_anno_group_redundance.xlsx",sep=""), overwrite = TRUE)
  }

}

RpsiAnnotator <- function(input_bedfile,
                          annotation_bed6file,
                          output_dir,
                          output_name)
{

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  psiAnnotator_path <- system.file("program", "psiAnnotator", package = "RpsiScan")
  system(paste("chmod 777", psiAnnotator_path, sep=" "))
  output_dir_output_name<-paste(output_dir,"/",output_name,sep="")

  input_bedfile6<-paste(output_dir,"/",basename(sub(".bed",".bed6",input_bedfile)),sep="")
  anno_bed<-paste(output_dir_output_name,"_anno.bed",sep="")
  anno_append_bed<-paste(output_dir_output_name,"_anno_append.bed",sep="")

  #cut -f 1-6 ${bedfile} > ${bedfile}6
  cmd0 <- paste("cut",
                "-f",
                "1-6",input_bedfile,
                ">", input_bedfile6,
                sep=" ")
  print(cmd0)
  system(cmd0)
  assign(paste(output_name,"_input_bed6_R",sep=""),fread(input_bedfile6),envir = .GlobalEnv)

  #bedAnnotator_cmd1="./script/bedAnnotator -s 1 --anno $bed6file --bed ${bedfile}6 -o ${bedfile%.bed}_anno.bed"
  #eval $bedAnnotator_cmd1 &
  cmd1 <- paste(psiAnnotator_path,
                "-s 1",
                "--anno", annotation_bed6file,
                "--bed", input_bedfile6,
                "-o", anno_bed,
                sep=" ")
  print(cmd1)
  system(cmd1)
  assign(paste(output_name,"_anno_bed_R",sep=""),fread(anno_bed,fill=TRUE),envir = .GlobalEnv)

  #bedAnnotator_cmd2="./script/bedAnnotator -s 1 --anno $bed6file --bed $bedfile -o ${bedfile%.bed}_anno_append.bed"
  #eval $bedAnnotator_cmd2 &
  cmd2 <- paste(psiAnnotator_path,
                "-s 1",
                "--anno", annotation_bed6file,
                "--bed", input_bedfile,
                "-o", anno_append_bed,
                sep=" ")
  print(cmd2)
  system(cmd2)
  assign(paste(output_name,"_anno_append_bed_R",sep=""),fread(anno_append_bed,fill=TRUE),envir = .GlobalEnv)

  anno.biotype.bed.file<-add_biotype(anno_bed_file = anno_bed,output_name)

  remove_redundancy(anno.biotype.bed.file=anno.biotype.bed.file, input.bed.file=input_bedfile,outfile=paste(output_dir,"/",basename(sub(".bed","",input_bedfile)),sep=""))

  print(paste("RpsiAnnotator result in",output_dir,sep=" "))
}
