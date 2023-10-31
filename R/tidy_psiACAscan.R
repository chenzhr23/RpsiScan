# tidy_psiACAscan
#
# This is a function to organize psiACAscan result into tidy data table

#' This is a function to organize psiACAscan result into tidy data table
#' @param infile The input file generated from RpsiACAscan function
#' @param appendfile The path to bed6.append file generated from RpsiACAscan function
#' @param orphaninfoFile The path to known snoRNA information file
#' @param psiU.SingleSites.bed psiU.SingleSites.bed pseudouridylation sites in bed format, default we recommend hg38.psiU.SingleSites.bed
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
tidy_psiACAscan <- function(infile,
                            appendfile,
                            orphaninfoFile,
                            psiU.SingleSites.bed,
                            output_dir,
                            output_name)
{

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  output_dir_output_name<-paste(output_dir,"/",output_name,sep="")

  #read in data
  pseudoU.sites<-read.delim(infile,header = F,sep=",")
  pseudoU.sites.append<-read.table(appendfile)
  pseudoU.sites.orphan<-read.csv(orphaninfoFile)
  colnames(pseudoU.sites.append)<-c("chro_id","stopRatioFC","annotation")

  #indexing
  ACA_index<-grep(">.*",pseudoU.sites$V1)
  ACA_seq_index<-ACA_index+1
  ACA_seq_sturcuture_index<-ACA_index+2

  #get by index
  ACA_id_index<-as.data.frame(pseudoU.sites$V1[ACA_index])
  colnames(ACA_id_index)<-"ACA_id_index"
  ACA_id_index$ACA_id_index<-str_replace_all(ACA_id_index$ACA_id_index,">","")
  ACA_seq<-as.data.frame(pseudoU.sites$V1[ACA_seq_index])
  colnames(ACA_seq)<-"ACA_seq"
  ACA_seq_sturcuture<-as.data.frame(pseudoU.sites$V1[ACA_seq_sturcuture_index])
  colnames(ACA_seq_sturcuture)<-"ACA_seq_sturcuture"

  search_df <- data.frame(ACA_id=ACA_id_index$ACA_id_index ,ACA_seq= ACA_seq$ACA_seq, ACA_seq_sturcuture= ACA_seq_sturcuture$ACA_seq_sturcuture )

  #ACA_id
  ACA_id<-str_extract_all(pseudoU.sites$V1,">.*")
  invisible(
    lapply(seq_along(ACA_id),function(x){
      if (length(ACA_id[[x]]) == 0) {
        ACA_id[[x]]<<-NA
      }
    }))
  ACA_id<-as.data.frame(unlist(ACA_id))
  colnames(ACA_id)<-"ACA_id"

  #target_site
  target_site<-str_extract_all(pseudoU.sites$V1,"#Target Site : .*")
  invisible(
    lapply(seq_along(target_site),function(x){
      if (length(target_site[[x]]) == 0) {
        target_site[[x]]<<-NA
      }
    }))
  # target_site<-target_site[-1]
  target_site<-as.data.frame(unlist(target_site))
  colnames(target_site)<-"target_site"
  # target_site[length(target_site$target_site)+1,]<-NA

  #target_score
  score_pairnum<-str_extract_all(pseudoU.sites$V1," Target Score: .*")
  invisible(lapply(seq_along(score_pairnum),function(x){
    if (length(score_pairnum[[x]]) == 0) {
      score_pairnum[[x]]<<-NA
    }
  }))

  score_pairnum<-as.data.frame(unlist(score_pairnum))
  colnames(score_pairnum)<-"score_pairnum"

  #upstream_sequence
  upstream_sequence<-str_extract_all(pseudoU.sites$V1," Upstream Sequence : .*")
  invisible(lapply(seq_along(upstream_sequence),function(x){
    if (length(upstream_sequence[[x]]) == 0) {
      upstream_sequence[[x]]<<-NA
    }
  }))
  upstream_sequence<-as.data.frame(unlist(upstream_sequence))
  colnames(upstream_sequence)<-"upstream_sequence"

  #upstream_structure
  upstream_structure<-str_extract_all(pseudoU.sites$V1," Upstream Structure: .*")
  invisible(lapply(seq_along(upstream_structure),function(x){
    if (length(upstream_structure[[x]]) == 0) {
      upstream_structure[[x]]<<-NA
    }
  }))
  upstream_structure<-as.data.frame(unlist(upstream_structure))
  colnames(upstream_structure)<-"upstream_structure"

  #target_sequence
  target_sequence<-str_extract_all(pseudoU.sites$V1," Target Sequence: .*")
  invisible(lapply(seq_along(target_sequence),function(x){
    if (length(target_sequence[[x]]) == 0) {
      target_sequence[[x]]<<-NA
    }
  }))
  target_sequence<-as.data.frame(unlist(target_sequence))
  colnames(target_sequence)<-"target_sequence"

  #forward
  forward<-str_extract_all(pseudoU.sites$V1," 5'.*3'")
  invisible(lapply(seq_along(forward),function(x){
    if (length(forward[[x]]) == 0) {
      # Do something
      forward[[x]]<<-NA
    }
  }))
  forward<-as.data.frame(unlist(forward))
  colnames(forward)<-"forward"

  #complement
  complement<-str_extract_all(pseudoU.sites$V1,".*\\|.*")
  invisible(lapply(seq_along(complement),function(x){
    if (length(complement[[x]]) == 0) {
      # Do something
      complement[[x]]<<-NA
    }
  }))
  complement<-as.data.frame(unlist(complement))
  colnames(complement)<-"complement"

  #reverse
  reverse<-str_extract_all(pseudoU.sites$V1," 3'.*5'")
  invisible(lapply(seq_along(reverse),function(x){
    if (length(reverse[[x]]) == 0) {
      # Do something
      reverse[[x]]<<-NA
    }
  }))
  reverse<-as.data.frame(unlist(reverse))
  colnames(reverse)<-"reverse"


  ACA_id<-ACA_id$ACA_id
  ACA_id<-str_replace_all(ACA_id,">","")
  pseudoU<-data.frame(ACA_id=ACA_id)

  target_site$target_site<-str_replace_all(target_site$target_site,"#Target Site : ","")
  Target_score_Modified_score_PairNum<-str_split(score_pairnum$score_pairnum,"\t")
  Target_score<-do.call(rbind.data.frame, lapply(Target_score_Modified_score_PairNum, `[`, 1))
  colnames(Target_score)<-"Target_score"
  Modified_score<-do.call(rbind.data.frame, lapply(Target_score_Modified_score_PairNum, `[`, 2))
  colnames(Modified_score)<-"Modified_score"
  PairNum<-do.call(rbind.data.frame, lapply(Target_score_Modified_score_PairNum, `[`, 3))
  colnames(PairNum)<-"PairNum"

  Target_score$Target_score<-str_replace_all(Target_score$Target_score," Target Score: ","")
  Modified_score$Modified_score<-str_replace_all(Modified_score$Modified_score,"Modified Score: ","")
  PairNum$PairNum<-str_replace_all(PairNum$PairNum,"PairNum: ","")

  upstream_sequence$upstream_sequence<-str_replace_all(upstream_sequence$upstream_sequence," Upstream Sequence : ","")
  upstream_structure$upstream_structure<-str_replace_all(upstream_structure$upstream_structure," Upstream Structure: ","")
  upstream_structure$upstream_structure<-str_replace_all(upstream_structure$upstream_structure,"\t"," ")
  target_sequence$target_sequence<-str_replace_all(target_sequence$target_sequence," Target Sequence: ","")

  pseudoU$target_site<-target_site$target_site
  pseudoU$Target_score<-Target_score$Target_score
  pseudoU$Modified_score<-Modified_score$Modified_score
  pseudoU$PairNum<-PairNum$PairNum
  pseudoU$upstream_sequence<-upstream_sequence$upstream_sequence
  pseudoU$upstream_structure<-upstream_structure$upstream_structure
  pseudoU$target_sequence<-target_sequence$target_sequence
  pseudoU$forward<-forward$forward
  pseudoU$complement<-complement$complement
  pseudoU$reverse<-reverse$reverse

  pseudoU<-as.data.frame(pseudoU)
  #edit(pseudoU)
  pseudoU<-pseudoU %>% fill(ACA_id)
  pseudoU<-pseudoU %>% mutate_at(c("Target_score"), funs(lead), n = 1 )
  pseudoU<-pseudoU %>% mutate_at(c("Modified_score"), funs(lead), n = 1 )
  pseudoU<-pseudoU %>% mutate_at(c("PairNum"), funs(lead), n = 1 )
  pseudoU<-pseudoU %>% mutate_at(c("upstream_sequence"), funs(lead), n = 2 )
  pseudoU<-pseudoU %>% mutate_at(c("upstream_structure"), funs(lead), n = 3 )
  pseudoU<-pseudoU %>% mutate_at(c("target_sequence"), funs(lead), n = 5 )
  pseudoU<-pseudoU %>% mutate_at(c("forward"), funs(lead), n = 6 )
  pseudoU<-pseudoU %>% mutate_at(c("complement"), funs(lead), n = 7 )
  pseudoU<-pseudoU %>% mutate_at(c("reverse"), funs(lead), n = 8 )

  pseudoU<-pseudoU[!is.na(pseudoU$target_site),]
  pseudoU<-pseudoU[!duplicated(pseudoU), ]
  pseudoU<-arrange(pseudoU,desc(Target_score))
  pseudoU<-pseudoU %>% separate(target_site,c("chro_id","target"),sep="\t")
  pseudoU<-pseudoU %>% left_join(pseudoU.sites.append,by=c("chro_id"="chro_id"))
  pseudoU$target<-str_replace_all(pseudoU$target," ","")
  pseudoU<-pseudoU %>% left_join(search_df,by=c("ACA_id"="ACA_id"))
  pseudoU<-pseudoU %>% left_join(pseudoU.sites.orphan,by=c("ACA_id"="ACA_id"))
  pseudoU$Target_score<-as.numeric(pseudoU$Target_score)
  pseudoU$diffPair<-sapply( strsplit(pseudoU$upstream_structure," ",fixed=T), "[", 3 )
  pseudoU$leftUnPair<-sapply( strsplit(pseudoU$upstream_structure," ",fixed=T), "[", 5 )
  pseudoU$rightUnPair<-sapply( strsplit(pseudoU$upstream_structure," ",fixed=T), "[", 7 )
  pseudoU$total_UnPair<-as.numeric(pseudoU$diffPair)+as.numeric(pseudoU$leftUnPair)+as.numeric(pseudoU$rightUnPair)
  tmp<-strsplit(pseudoU$annotation,"|",fixed = T)
  pseudoU$enstid<-NULL
  pseudoU$enstsymbol<-NULL
  pseudoU$ensgid<-NULL
  pseudoU$ensgsymbol<-NULL
  pseudoU$type<-NULL
  invisible(
    lapply(seq_along(tmp), function(x){
      pseudoU$enstid[x]<<-tmp[[x]][1]
      pseudoU$enstsymbol[x]<<-tmp[[x]][2]
      pseudoU$ensgid[x]<<-tmp[[x]][3]
      pseudoU$ensgsymbol[x]<<-tmp[[x]][4]
      pseudoU$type[x]<<-tmp[[x]][5]
    }))
  target<-read.table(psiU.SingleSites.bed)
  target$uniq_id<-paste(target$V1,target$V2,target$V3,target$V6,sep="_")
  pseudoU<-pseudoU %>% left_join(target,by=c("chro_id"="uniq_id"))
  pseudoU<-subset(pseudoU, select=-c(V1,V2,V3,V5,V6))

  if(is.character(pseudoU[is.na(pseudoU$V4),]$V4)){
    print("all of the snoRNAs are with known snoRNA targets!")
  }else{
    pseudoU[is.na(pseudoU$V4),]$V4<-"unknown"
  }

  pseudoU<-pseudoU%>%rename(annotated_snoTar=V4)
  pseudoU$pair_id<-paste(pseudoU$ACA_id,pseudoU$chro_id,sep="_")
  orphan_na_index<-which(is.na(pseudoU$Orphan))
  pseudoU[orphan_na_index,]$Orphan<-"TRUE"

  orphan_FALSE_index<-which(pseudoU$annotated_snoTar!="unknown")
  if(is.integer(orphan_FALSE_index)){
    print("none of the snoRNAs are candidate orphan snoRNA!")
  }else{
    pseudoU[orphan_FALSE_index,]$Orphan<-"FALSE"
  }

  write.csv(pseudoU,paste(output_dir_output_name,"_pseudoU.csv",sep=""),row.names = F,quote = F)
  write.xlsx(pseudoU,paste(output_dir_output_name,"_pseudoU.xlsx",sep=""),overwrite=T)


  #pseudoU (overall) ##stopration/target_score bubble plot
  mycol=brewer.pal(9, "Set1")
  to_bubble<-pseudoU
  to_bubble <- to_bubble %>% mutate(
    isorphan = case_when(
      target_summary=="orphan" ~ "Y",
      target_summary!="orphan" ~ "N",
      TRUE ~  "U"
    )
  )
  to_bubble_select<-to_bubble %>% select(ACA_id,enstid,enstsymbol,Target_score,stopRatioFC,isorphan)
  to_bubble_select<-arrange(to_bubble_select,desc(Target_score))

  # if(length(to_bubble_select$ACA_id)>50){
  to_bubble_select<-head(to_bubble_select,200)
  if(dim(pseudoU)[1]>0){
    p_bubble<-ggplot(to_bubble_select, aes(x=ACA_id, y=paste(enstid,enstsymbol,sep="_"), size =stopRatioFC ,colour =Target_score,label =isorphan,fill=isorphan)) +
      ylab("annotation")+
      xlab("ACA_id")+
      geom_point()+
      geom_text(colour="black")+
      theme_bw()+
      theme(plot.margin = unit(c(2,2,2,2), "cm"))+
      # scale_color_gradient(low = mycol[2],high = mycol[1])+
      scale_color_gradient(low =brewer.pal(9, "Set3")[2] ,high = mycol[1])+
      theme(axis.text.x=element_text(angle=90,vjust=0.5))+scale_fill_manual(values=c("black","black","black"),labels=c("N-No","Y-Yes","U-Unknown"))

    p_bubble
    ggsave(paste(output_dir_output_name,"_pseudoU.pdf",sep=""),width = 20, height = 15)
    # }else{
    #   print("ACA target result less than 400, no bubble plot")
    # }
  }

  #pie chart for orphan snoRNA
  if(dim(pseudoU)[1]>0){
    pseudoU_Orphan<-arrange(as.data.frame(table(pseudoU$Orphan)),desc(Freq))
    percent<-round(100*pseudoU_Orphan$Freq/sum(pseudoU_Orphan$Freq),2)
    percent <-paste('(',percent, "%", ", ", pseudoU_Orphan$Freq,')', sep = "")
    pseudoU_Orphan$Var1 <- paste(pseudoU_Orphan$Var1, percent, sep = '')
    pseudoU_Orphan$Var1 <- factor(pseudoU_Orphan$Var1, levels = pseudoU_Orphan$Var1)
    mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"),brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Pastel1"),brewer.pal(8, "Pastel2"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))
    pie_plot<-ggplot(data = pseudoU_Orphan, mapping = aes(x = 'Content', y = Freq, fill = Var1)) +
      geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +
      theme(axis.text = element_blank()) +
      theme(axis.ticks = element_blank()) +
      theme(panel.grid = element_blank(),
            panel.background=element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank()) +
      guides(fill = guide_legend(title = "Orphan snoRNA-guided RNA(snoDB)"))+
      scale_fill_manual(values = mycol)

    pdf(paste(output_dir_output_name,"_pseudoU_Orphan_piechart.pdf",sep=""))
    print(pie_plot)
    invisible(dev.off())
  }

  #get high confidence result
  pseudoU$Modified_score<-as.numeric(pseudoU$Modified_score)
  pseudoU$PairNum<-as.numeric(pseudoU$PairNum)
  pseudoU_high_confidence<-pseudoU %>% filter(Modified_score>=15 & PairNum>=9 & total_UnPair<=2)#Modified Score:12; PairNum: 9; total_UnPair: <=2
  write.csv(pseudoU_high_confidence,paste(output_dir_output_name,"_pseudoU_high_confidence.csv",sep=""),row.names = F,quote = F)
  write.xlsx(pseudoU_high_confidence,paste(output_dir_output_name,"_pseudoU_high_confidence.xlsx",sep=""),overwrite=T)
  write.table(pseudoU_high_confidence,paste(output_dir_output_name,"_pseudoU_high_confidence.txt",sep=""),row.names = F,quote = F)

  #pie chart for all snoRNA-guided RNAs
  mRNA<-data.frame(gene_biotype="mRNA",feature_type=c("protein_coding"))
  pseudogene<-data.frame(gene_biotype="pseudogene",feature_type=c("rRNA_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_processed_pseudogene", "unitary_pseudogene", "pseudogene", "polymorphic_pseudogene", "transcribed_unitary_pseudogene", "TR_V_pseudogene", "TR_J_pseudogene", "IG_V_pseudogene", "IG_C_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "IG_J_pseudogene"))
  lncRNA<-data.frame(gene_biotype="lncRNA",feature_type=c("lncRNA","processed_transcript","lincRNA","non_coding","3prime_overlapping_ncRNA","3prime_overlapping_ncrna","sense_intronic","antisense","sense_overlapping","known_ncrna","macro_lncRNA","bidirectional_promoter_lncRNA","retained_intron","TEC"))
  sncRNA<-data.frame(gene_biotype="sncRNA",feature_type=c("snRNA","snoRNA","misc_RNA","miscRNA","miRNA","ribozyme","sRNA","scRNA","scaRNA","srpRNA","tRNA-Deu","tRNA-RTE","piRNA","siRNA"))
  rRNA<-data.frame(gene_biotype="rRNA",feature_type=c("rRNA","Mt_rRNA"))
  tRNA<-data.frame(gene_biotype="tRNA",feature_type=c("tRNA","Mt_tRNA","vaultRNA"))
  IG_gene<-data.frame(gene_biotype="IG_gene",feature_type=c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene"))
  TR_gene<-data.frame(gene_biotype="TR_gene",feature_type=c("TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene"))
  repeatMasker<-data.frame(gene_biotype="repeatMasker",feature_type=c("5S-Deu-L2","Alu","centr","CR1","DNA","DNA?","ERV1","ERV1?","ERVK","ERVL","ERVL?","ERVL-MaLR","Gypsy","Gypsy?","hAT","hAT?","hAT-Ac","hAT-Blackjack","hAT-Charlie","hAT-Tag1","hAT-Tip100","hAT-Tip100?","Helitron","Helitron?","L1","L2","Low_complexity","LTR","LTR?","MIR","MULE-MuDR","nonsense_mediated_decay","non_stop_decay","Penelope","PiggyBac","PiggyBac?","RNA","RTE-BovB","RTE-X","Satellite","Simple_repeat","SVA","TcMar?","TcMar-Mariner","TcMar-Tc2","TcMar-Tigger","telo","Unknown","acro","Crypton","Dong-R4","I-Jockey","Kolobok","L1-Tx1","Merlin","MULE-MuDR?","PIF-Harbinger","SINE?","TcMar","TcMar-Pogo","TcMar-Tc1"))
  intergenic<-data.frame(gene_biotype="intergenic",feature_type=c("intergenic"))
  circRNA<-data.frame(gene_biotype="circRNA",feature_type=c("circRNA"))
  category<-rbind(mRNA,pseudogene,lncRNA,sncRNA,rRNA,tRNA,IG_gene,TR_gene,repeatMasker,intergenic,circRNA)

  # pseudoU_uniq_chro_id<-pseudoU_high_confidence %>% distinct(chro_id,.keep_all = TRUE)
  pseudoU_uniq_chro_id<-pseudoU_high_confidence %>% distinct(target_sequence,.keep_all = TRUE)
  pseudoU_uniq_chro_id_add_gene_biotype<-pseudoU_uniq_chro_id %>% left_join(category,by=c("type"="feature_type"))

  pseudoU_RNA_type<-arrange(as.data.frame(table(pseudoU_uniq_chro_id_add_gene_biotype$gene_biotype)),desc(Freq))
  percent<-round(100*pseudoU_RNA_type$Freq/sum(pseudoU_RNA_type$Freq),2)
  percent <-paste('(',percent, "%", ", ", pseudoU_RNA_type$Freq,')', sep = "")
  pseudoU_RNA_type$Var1 <- paste(pseudoU_RNA_type$Var1, percent, sep = '')
  pseudoU_RNA_type$Var1 <- factor(pseudoU_RNA_type$Var1, levels = pseudoU_RNA_type$Var1)
  mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"),brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Pastel1"),brewer.pal(8, "Pastel2"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))
  pie_plot<-ggplot(data = pseudoU_RNA_type, mapping = aes(x = 'Content', y = Freq, fill = Var1)) +
    geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(panel.grid = element_blank(),
          panel.background=element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()) +
    guides(fill = guide_legend(title = "snoRNA-guided RNA"))+
    scale_fill_manual(values = mycol)

  pdf(paste(output_dir_output_name,"_pseudoU_high_confidence_snoRNA_guided_RNA_piechart.pdf",sep=""))
  print(pie_plot)
  invisible(dev.off())

  pseudoU$snoRNA_target_duplex<-paste("\n",pseudoU$forward, "\n",pseudoU$complement,"\n", pseudoU$reverse,"\n",sep="")
  pseudoU<-subset(pseudoU, select=-c(forward,complement,reverse))
  write.table(pseudoU,paste(output_dir_output_name,"_pseudoU.txt",sep=""),row.names = F,quote = F)

  print(paste("tidy_psiACAscan result in",output_dir,sep=" "))
}
