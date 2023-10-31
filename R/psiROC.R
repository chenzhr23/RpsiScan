# psiROC
#
# Perform ROC evaluation based on ground truth pseudouridylation data set

#' Perform ROC evaluation based on ground truth pseudouridylation data set
#' @param rocfile ROC input file of single sites information (with suffix _roc_plot.txt) generated from ground_truth() function
#' @param rRNAfile default we recommend hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
#' @param rRNAfile2 default we recommend hg38.psiU.SingleSites.bed
#' @param filtfile File to be filted
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
#'

#font function from ggpubr
font<-function (object, size = NULL, color = NULL, face = NULL, family = NULL,
                ...)
{
  elmt <- element_text(size = size, color = color, face = face,
                       family = family, ...)
  switch(object, title = theme(plot.title = elmt), subtitle = theme(plot.subtitle = elmt),
         caption = theme(plot.caption = elmt), x = theme(axis.title.x = elmt),
         xlab = theme(axis.title.x = elmt), x.title = theme(axis.title.x = elmt),
         y = theme(axis.title.y = elmt), ylab = theme(axis.title.y = elmt),
         y.title = theme(axis.title.y = elmt), xy = theme(axis.title.x = elmt,
                                                          axis.title.y = elmt), xylab = theme(axis.title.x = elmt,
                                                                                              axis.title.y = elmt), xy.title = theme(axis.title.x = elmt,
                                                                                                                                     axis.title.y = elmt), axis.title = theme(axis.title.x = elmt,
                                                                                                                                                                              axis.title.y = elmt), legendtitle = theme(legend.title = elmt),
         legend.title = theme(legend.title = elmt), legendtext = theme(legend.text = elmt),
         legend.text = theme(legend.text = elmt), x.text = theme(axis.text.x = elmt),
         y.text = theme(axis.text.y = elmt), xy.text = theme(axis.text.x = elmt,
                                                             axis.text.y = elmt), yxtext = theme(axis.text.x = elmt,
                                                                                                 axis.text.y = elmt), axis.text = theme(axis.text.x = elmt,
                                                                                                                                        axis.text.y = elmt), stop("Don't support ", object))
}

#.method_info function from ggpubr
.method_info<-function (method)
{
  if (is.null(method))
    method = "wilcox.test"
  allowed.methods <- list(t = "t.test", t.test = "t.test",
                          student = "t.test", wiloxon = "wilcox.test", wilcox = "wilcox.test",
                          wilcox.test = "wilcox.test", anova = "anova", aov = "anova",
                          kruskal = "kruskal.test", kruskal.test = "kruskal.test")
  method.names <- list(t.test = "T-test", wilcox.test = "Wilcoxon",
                       anova = "Anova", kruskal.test = "Kruskal-Wallis")
  if (!(method %in% names(allowed.methods)))
    stop("Non-supported method specified. Allowed methods are one of: ",
         .collapse(allowed.methods, sep = ", "))
  method <- allowed.methods[[method]]
  method.name <- method.names[[method]]
  list(method = method, name = method.name)
}

#.add_item function from ggpubr
.add_item<-function (.list, ...)
{
  pms <- list(...)
  for (pms.names in names(pms)) {
    .list[[pms.names]] <- pms[[pms.names]]
  }
  .list
}

#.is_p.signif_in_mapping function from ggpubr
.is_p.signif_in_mapping<-function (mapping)
{
  res <- FALSE
  if (!is.null(mapping)) {
    if (!is.null(mapping$label)) {
      .label <- rlang::as_label(mapping$label)
      res <- grepl(pattern = "p\\.signif", .label)
    }
  }
  return(res)
}

#.is_empty function from ggpubr
.is_empty<-function (x)
{
  length(x) == 0
}

#stat_compare_means function from ggpubr
stat_compare_means<-function (mapping = NULL, data = NULL, method = NULL, paired = FALSE,
                              method.args = list(), ref.group = NULL, comparisons = NULL,
                              hide.ns = FALSE, label.sep = ", ", label = NULL, label.x.npc = "left",
                              label.y.npc = "top", label.x = NULL, label.y = NULL, vjust = 0,
                              tip.length = 0.03, bracket.size = 0.3, step.increase = 0,
                              symnum.args = list(), geom = "text", position = "identity",
                              na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...)
{
  if (!is.null(comparisons)) {
    method.info <- .method_info(method)
    method <- method.info$method
    method.args <- .add_item(method.args, paired = paired)
    pms <- list(...)
    size <- ifelse(is.null(pms$size), 3.88, pms$size)
    color <- ifelse(is.null(pms$color), "black", pms$color)
    map_signif_level <- FALSE
    if (is.null(label))
      label <- "p.format"
    if (.is_p.signif_in_mapping(mapping) | (label %in% "p.signif")) {
      map_signif_level <- c(`****` = 1e-04, `***` = 0.001,
                            `**` = 0.01, `*` = 0.05, ns = Inf)
      if (hide.ns)
        map_signif_level <- .hide_ns(map_signif_level)
    }
    if (!.is_empty(symnum.args)) {
      symnum.args.isok <- length(symnum.args$cutpoints ==
                                   length(symnum.args$symbols))
      if (!symnum.args.isok)
        stop("Incorrect format detected in symnum.args. ",
             "Check the documentation.")
      map_signif_level <- symnum.args$cutpoints[-1]
      names(map_signif_level) <- symnum.args$symbols
      if (hide.ns)
        map_signif_level <- .hide_ns(map_signif_level)
    }
    if (missing(step.increase)) {
      step.increase <- ifelse(is.null(label.y), 0.12, 0)
    }
    ggsignif::geom_signif(comparisons = comparisons, y_position = label.y,
                          test = method, test.args = method.args, step_increase = step.increase,
                          size = bracket.size, textsize = size, color = color,
                          map_signif_level = map_signif_level, tip_length = tip.length,
                          data = data, vjust = vjust)
  }
  else {
    mapping <- .update_mapping(mapping, label)
    layer(stat = StatCompareMeans, data = data, mapping = mapping,
          geom = geom, position = position, show.legend = show.legend,
          inherit.aes = inherit.aes, params = list(label.x.npc = label.x.npc,
                                                   label.y.npc = label.y.npc, label.x = label.x,
                                                   label.y = label.y, label.sep = label.sep, method = method,
                                                   method.args = method.args, paired = paired, ref.group = ref.group,
                                                   symnum.args = symnum.args, hide.ns = hide.ns,
                                                   na.rm = na.rm, vjust = vjust, ...))
  }
}

psiROC <- function(rocfile, rRNAfile, rRNAfile2, filtfile,output_dir,output_name)
{
  #nohup Rscript $(dirname "$0")/roc.r -f ${outFileName}_roc_plot.txt -t ${outFileName}_roc_filt.bed -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed -o ${outFileName} > ${outFileName}_roc_bestthres.log 2>&1 &

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  outfile_prefix<-paste(output_dir,"/",output_name,sep="")

  # load roc_plot.txt
  ROC_data<-as.data.frame(fread(rocfile))
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
                        "class")
  # ROC_data$base<-str_replace_all(as.character(ROC_data$base),"TRUE","T")
  ROC_data_sel<-ROC_data %>% select(tretRtsRatio,rtsRatioFold,tretBefRatio,befRatioFold,tretAftRatio,aftRatioFold,tretMutRatio,mutRatioFold,tretDelRatio,delRatioFold,class)


  #rRNA-psi-non-psi visualization
  ROC_data_melt<-melt(ROC_data_sel,id.vars = "class")
  ROC_data_melt$class<-str_replace(as.character(ROC_data_melt$class),"0","non-psi")
  ROC_data_melt$class<-str_replace(as.character(ROC_data_melt$class),"1","psi")

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


  cat("=====================tretRtsRatio/ctrlRtsRatio=========================\n")
  #tretRtsRatio ctrlRtsRatio
  input=roc(ROC_data$class,ROC_data$ctrlRtsRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # pdf(paste(outfile_prefix, '_roc_tretRtsRatio.pdf', sep=""))
  CMC=roc(ROC_data$class,ROC_data$tretRtsRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # invisible(dev.off())
  input_label = paste("ctrlRtsRatio AUC:",sprintf("%.3f",input$auc))
  CMC_label = paste("tretRtsRatio AUC:",sprintf("%.3f",CMC$auc))
  preList=list(CMC=CMC, Input=input)
  names(preList) <- c(CMC_label,input_label)
  RtsRatio_plot<-ggroc(preList,size=0.8)
  RtsRatio_plot<-RtsRatio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  RtsRatio_plot<-RtsRatio_plot+ guides(colour = guide_legend(nrow = 2))
  RtsRatio_plot<-RtsRatio_plot+ annotate("text", x = .5, y = .5,
                                         label = paste("Method: ",roc.test(input, CMC)$method,"\np.value: ",signif(roc.test(input, CMC)$p.value,3),sep=""),
                                         size = 3)
  print(paste(outfile_prefix,"_RtsRatio_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_RtsRatio_plot.pdf",sep=""),width =4 ,height = 3)
  RtsRatio_plot
  invisible(dev.off())
  roc.test(input, CMC)


  cat("=====================tretBefRatio/ctrlBefRatio=========================\n")
  #tretBefRatio controlAfterRpmFold
  input=roc(ROC_data$class,ROC_data$ctrlBefRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # pdf(paste(outfile_prefix, '_roc_tretBefRatio.pdf', sep=""))
  CMC=roc(ROC_data$class,ROC_data$tretBefRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # invisible(dev.off())
  input_label = paste("ctrlBefRatio AUC:",sprintf("%.3f",input$auc))
  CMC_label = paste("tretBefRatio AUC:",sprintf("%.3f",CMC$auc))
  aftList=list(CMC=CMC, Input=input)
  names(aftList) <- c(CMC_label,input_label)
  BefRatio_plot<-ggroc(aftList,size=0.8)
  BefRatio_plot<-BefRatio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  BefRatio_plot<-BefRatio_plot+ guides(colour = guide_legend(nrow = 2))
  BefRatio_plot<-BefRatio_plot+ annotate("text", x = .5, y = .5,
                                         label = paste("Method: ",roc.test(input, CMC)$method,"\np.value: ",signif(roc.test(input, CMC)$p.value,3),sep=""),
                                         size = 3)
  print(paste(outfile_prefix,"_BefRatio_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_BefRatio_plot.pdf",sep=""),width =4 ,height = 3)
  BefRatio_plot
  invisible(dev.off())
  roc.test(input, CMC)


  cat("=====================tretAftRatio/ctrlAftRatio=========================\n")
  #tretAftRatio ctrlAftRatio
  input=roc(ROC_data$class,ROC_data$ctrlAftRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # pdf(paste(outfile_prefix, '_roc_tretAftRatio.pdf', sep=""))
  CMC=roc(ROC_data$class,ROC_data$tretAftRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # invisible(dev.off())
  ctrlAftRatio_label = paste("ctrlAftRatio AUC:",sprintf("%.3f",input$auc))
  tretAftRatio_label = paste("tretAftRatio AUC:",sprintf("%.3f",CMC$auc))
  StopMeanFoldList=list(CMC=CMC, Input=input)
  names(StopMeanFoldList) <- c(tretAftRatio_label,ctrlAftRatio_label)
  AftRatio_plot<-ggroc(StopMeanFoldList,size=0.8)
  AftRatio_plot<-AftRatio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  AftRatio_plot<-AftRatio_plot+ guides(colour = guide_legend(nrow = 2))
  AftRatio_plot<-AftRatio_plot+ annotate("text", x = .5, y = .5,
                                         label = paste("Method: ",roc.test(input, CMC)$method,"\np.value: ",signif(roc.test(input, CMC)$p.value,3),sep=""),
                                         size = 3)
  print(paste(outfile_prefix,"_AftRatio_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_AftRatio_plot.pdf",sep=""),width =4 ,height = 3)
  AftRatio_plot
  invisible(dev.off())

  cat("=====================tretMutRatio/ctrlMutRatio=========================\n")
  #tretMutRatio ctrlMutRatio
  input=roc(ROC_data$class,ROC_data$ctrlMutRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # pdf(paste(outfile_prefix, '_roc_tretMutRatio.pdf', sep=""))
  CMC=roc(ROC_data$class,ROC_data$tretMutRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # invisible(dev.off())
  ctrlMutRatio_label = paste("ctrlMutRatio AUC:",sprintf("%.3f",input$auc))
  tretMutRatio_label = paste("tretMutRatio AUC:",sprintf("%.3f",CMC$auc))
  StopMeanFoldList=list(CMC=CMC, Input=input)
  names(StopMeanFoldList) <- c(tretMutRatio_label,ctrlMutRatio_label)
  MutRatio_plot<-ggroc(StopMeanFoldList,size=0.8)
  MutRatio_plot<-MutRatio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  MutRatio_plot<-MutRatio_plot+ guides(colour = guide_legend(nrow = 2))
  MutRatio_plot<-MutRatio_plot+ annotate("text", x = .5, y = .5,
                                         label = paste("Method: ",roc.test(input, CMC)$method,"\np.value: ",signif(roc.test(input, CMC)$p.value,3),sep=""),
                                         size = 3)
  print(paste(outfile_prefix,"_MutRatio_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_MutRatio_plot.pdf",sep=""),width =4 ,height = 3)
  MutRatio_plot
  invisible(dev.off())
  roc.test(input, CMC)

  cat("=====================tretDelRatio/ctrlDelRatio=========================\n")
  #tretDelRatio ctrlDelRatio
  input=roc(ROC_data$class,ROC_data$ctrlDelRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # pdf(paste(outfile_prefix, '_roc_tretDelRatio.pdf', sep=""))
  CMC=roc(ROC_data$class,ROC_data$tretDelRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  # invisible(dev.off())
  ctrlDelRatio_label = paste("ctrlDelRatio AUC:",sprintf("%.3f",input$auc))
  tretDelRatio_label = paste("tretDelRatio AUC:",sprintf("%.3f",CMC$auc))
  StopMeanFoldList=list(CMC=CMC, Input=input)
  names(StopMeanFoldList) <- c(tretDelRatio_label,ctrlDelRatio_label)
  DelRatio_plot<-ggroc(StopMeanFoldList,size=0.8)
  DelRatio_plot<-DelRatio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  DelRatio_plot<-DelRatio_plot+ guides(colour = guide_legend(nrow = 2))
  DelRatio_plot<-DelRatio_plot+ annotate("text", x = .5, y = .5,
                                         label = paste("Method: ",roc.test(input, CMC)$method,"\np.value: ",signif(roc.test(input, CMC)$p.value,3),sep=""),
                                         size = 3)
  print(paste(outfile_prefix,"_DelRatio_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_DelRatio_plot.pdf",sep=""),width =4 ,height = 3)
  DelRatio_plot
  invisible(dev.off())

  cat("=====================ten variables: tretRtsRatio, rtsRatioFold, tretBefRatio, befRatioFold, tretAftRatio, aftRatioFold,tretMutRatio, mutRatioFold, tretDelRatio, delRatioFold=========================\n")
  mycol = brewer.pal(10, "Paired")[c(1:10)]
  pdf(paste(outfile_prefix, '_roc_tretRtsRatio.pdf', sep=""))
  tretRtsRatio=roc(ROC_data$class,ROC_data$tretRtsRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_rtsRatioFold.pdf', sep=""))
  rtsRatioFold=roc(ROC_data$class,ROC_data$rtsRatioFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_tretBefRatio.pdf', sep=""))
  tretBefRatio=roc(ROC_data$class,ROC_data$tretBefRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_befRatioFold.pdf', sep=""))
  befRatioFold=roc(ROC_data$class,ROC_data$befRatioFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_tretAftRatio.pdf', sep=""))
  tretAftRatio=roc(ROC_data$class,ROC_data$tretAftRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_aftRatioFold.pdf', sep=""))
  aftRatioFold=roc(ROC_data$class,ROC_data$aftRatioFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_tretMutRatio.pdf', sep=""))
  tretMutRatio=roc(ROC_data$class,ROC_data$tretMutRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_mutRatioFold.pdf', sep=""))
  mutRatioFold=roc(ROC_data$class,ROC_data$mutRatioFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_tretDelRatio.pdf', sep=""))
  tretDelRatio=roc(ROC_data$class,ROC_data$tretDelRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_delRatioFold.pdf', sep=""))
  delRatioFold=roc(ROC_data$class,ROC_data$delRatioFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  invisible(dev.off())
  tretRtsRatio_label = paste("tretRtsRatio AUC:",sprintf("%.3f",tretRtsRatio$auc))
  rtsRatioFold_label = paste("rtsRatioFold AUC:",sprintf("%.3f",rtsRatioFold$auc))
  tretBefRatio_label = paste("tretBefRatio AUC:",sprintf("%.3f",tretBefRatio$auc))
  befRatioFold_label = paste("befRatioFold AUC:",sprintf("%.3f",befRatioFold$auc))
  tretAftRatio_label = paste("tretAftRatio AUC:",sprintf("%.3f",tretAftRatio$auc))
  aftRatioFold_label = paste("aftRatioFold AUC:",sprintf("%.3f",aftRatioFold$auc))
  tretMutRatio_label = paste("tretMutRatio AUC:",sprintf("%.3f",tretMutRatio$auc))
  mutRatioFold_label = paste("mutRatioFold AUC:",sprintf("%.3f",mutRatioFold$auc))
  tretDelRatio_label = paste("tretDelRatio AUC:",sprintf("%.3f",tretDelRatio$auc))
  delRatioFold_label = paste("delRatioFold AUC:",sprintf("%.3f",delRatioFold$auc))
  ratioList=list(tretRtsRatio = tretRtsRatio, rtsRatioFold = rtsRatioFold, tretBefRatio = tretBefRatio, befRatioFold = befRatioFold, tretAftRatio = tretAftRatio, aftRatioFold= aftRatioFold, tretMutRatio = tretMutRatio,  mutRatioFold = mutRatioFold, tretDelRatio = tretDelRatio, delRatioFold = delRatioFold)
  names(ratioList) <- c(tretRtsRatio_label, rtsRatioFold_label, tretBefRatio_label, befRatioFold_label, tretAftRatio_label, aftRatioFold_label,tretMutRatio_label,mutRatioFold_label,tretDelRatio_label,delRatioFold_label)
  ten_variables_plot<-ggroc(ratioList,size=0.8)
  ten_variables_plot<-ten_variables_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  ten_variables_plot<-ten_variables_plot+ guides(colour = guide_legend(nrow = 5, ncol = 2))+scale_color_manual(values = mycol)
  print(paste(outfile_prefix,"_ten_variables_ratio_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_ten_variables_plot.pdf",sep=""),width =8 ,height = 7)
  ten_variables_plot
  invisible(dev.off())

  pdf(paste(outfile_prefix,"_roc_summary.pdf",sep=""),width=13,height=8)
  plot_grid(RtsRatio_plot, BefRatio_plot, AftRatio_plot, MutRatio_plot, DelRatio_plot,ten_variables_plot, align = "hv",labels = c('A','B','C','D','E','F'))
  invisible(dev.off())

  invisible(lapply(seq_along(ROC_data_sel[,-length(ROC_data_sel)]),function(x){
    assign(colnames(ROC_data_sel)[x],roc(ROC_data$class,ROC_data_sel[,x],smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T),envir = .GlobalEnv)
    assign(paste(colnames(ROC_data_sel)[x],"_thres",sep=""),coords(get(colnames(ROC_data_sel)[x]), "best",transpose = TRUE)[1],envir = .GlobalEnv)
  }))

  # cat("=====================tretRtsRatio threshold=========================\n")
  # tretRtsRatio_thres<-coords(CMC, "best",transpose = TRUE)[1]
  # coords(tretRtsRatio, "best", ret = "all", transpose = TRUE)


  thres<-cbind(tretRtsRatio_thres,rtsRatioFold_thres,tretBefRatio_thres,befRatioFold_thres,tretAftRatio_thres, aftRatioFold_thres, tretMutRatio_thres, mutRatioFold_thres, tretDelRatio_thres, delRatioFold_thres)
  thres
  write.table(thres,paste(outfile_prefix, '_roc_bestthres.txt', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)
  write.table(thres,paste(outfile_prefix, '_roc_bestthres_colname.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)

  # Values of all the confusion matrix terms were calculated at the optimal threshold
  invisible(lapply(seq_along(ROC_data_sel[,-length(ROC_data_sel)]),function(x){
    assign(paste(colnames(ROC_data_sel)[x],"_TP",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==1 & ROC_data_sel[,x] > get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#True Positives (TP)
    assign(paste(colnames(ROC_data_sel)[x],"_FP",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==0 & ROC_data_sel[,x] > get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#False Positives (FP)
    assign(paste(colnames(ROC_data_sel)[x],"_TN",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==0 & ROC_data_sel[,x] <= get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#True Negatives (TN)
    assign(paste(colnames(ROC_data_sel)[x],"_FN",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==1 & ROC_data_sel[,x] <= get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#False Negatives (FN)
    assign(paste(colnames(ROC_data_sel)[x],"_TPR",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FN",sep=""))),3),envir = .GlobalEnv)#sensitivity (true positive rate, TPR)
    assign(paste(colnames(ROC_data_sel)[x],"_TNR",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TN",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TN",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FP",sep=""))),3),envir = .GlobalEnv)#specifity (selectivity or true negative rate, TNR)
    assign(paste(colnames(ROC_data_sel)[x],"_FPR",sep=""),round(1-get(paste(colnames(ROC_data_sel)[x],"_TNR",sep="")),3),envir = .GlobalEnv)#False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
    assign(paste(colnames(ROC_data_sel)[x],"_FNR",sep=""),round(1-get(paste(colnames(ROC_data_sel)[x],"_TPR",sep="")),3),envir = .GlobalEnv)#False Negative Rate, FNR)
    assign(paste(colnames(ROC_data_sel)[x],"_Prec",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FP",sep=""))),3),envir = .GlobalEnv)#Precision
    assign(paste(colnames(ROC_data_sel)[x],"_Recall",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FN",sep=""))),3),envir = .GlobalEnv)#Recall
    assign(paste(colnames(ROC_data_sel)[x],"_ACC",sep=""),round((get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_TN",sep="")))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_TN",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FN",sep=""))),3),envir = .GlobalEnv)#accuracy
    assign(paste(colnames(ROC_data_sel)[x],"_F1_score",sep=""),round((2*get(paste(colnames(ROC_data_sel)[x],"_Recall",sep=""))*get(paste(colnames(ROC_data_sel)[x],"_Prec",sep="")))/(get(paste(colnames(ROC_data_sel)[x],"_Recall",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_Prec",sep=""))),3),envir = .GlobalEnv)#F1_score
    assign(paste(colnames(ROC_data_sel)[x],"_eval",sep=""),cbind(get(paste(colnames(ROC_data_sel)[x],"_TP",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FP",sep="")),get(paste(colnames(ROC_data_sel)[x],"_TN",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FN",sep="")),get(paste(colnames(ROC_data_sel)[x],"_TPR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_TNR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FPR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FNR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_Prec",sep="")),get(paste(colnames(ROC_data_sel)[x],"_Recall",sep="")),get(paste(colnames(ROC_data_sel)[x],"_ACC",sep="")),get(paste(colnames(ROC_data_sel)[x],"_F1_score",sep=""))),envir = .GlobalEnv)
    tmp<-as.data.frame(get(paste(colnames(ROC_data_sel)[x],"_eval",sep="")))
    names(tmp)<-c(
      paste(colnames(ROC_data_sel)[x],"_TP",sep=""),
      paste(colnames(ROC_data_sel)[x],"_FP",sep=""),
      paste(colnames(ROC_data_sel)[x],"_TN",sep=""),
      paste(colnames(ROC_data_sel)[x],"_FN",sep=""),
      paste(colnames(ROC_data_sel)[x],"_TPR",sep=""),
      paste(colnames(ROC_data_sel)[x],"_TNR",sep=""),
      paste(colnames(ROC_data_sel)[x],"_FPR",sep=""),
      paste(colnames(ROC_data_sel)[x],"_FNR",sep=""),
      paste(colnames(ROC_data_sel)[x],"_Prec",sep=""),
      paste(colnames(ROC_data_sel)[x],"_Recall",sep=""),
      paste(colnames(ROC_data_sel)[x],"_ACC",sep=""),
      paste(colnames(ROC_data_sel)[x],"_F1_score",sep="")
    )
    assign(paste(colnames(ROC_data_sel)[x],"_eval",sep=""), tmp, envir = .GlobalEnv)
  }))


  confusion_matrix_and_indicators<-as.data.frame(rbind(
    t(as.data.frame(tretRtsRatio_eval)),
    t(as.data.frame(rtsRatioFold_eval)),
    t(as.data.frame(tretBefRatio_eval)),
    t(as.data.frame(befRatioFold_eval)),
    t(as.data.frame(tretAftRatio_eval)),
    t(as.data.frame(aftRatioFold_eval)),
    t(as.data.frame(tretMutRatio_eval)),
    t(as.data.frame(mutRatioFold_eval)),
    t(as.data.frame(tretDelRatio_eval)),
    t(as.data.frame(delRatioFold_eval))))
  colnames(confusion_matrix_and_indicators)<-"value_or_percentage"
  confusion_matrix_and_indicators_arrange<-arrange(confusion_matrix_and_indicators,desc(value_or_percentage))
  write.table(confusion_matrix_and_indicators,paste(outfile_prefix, '_roc_confusion_matrix_and_indicators.txt', sep=""),sep="\t",row.names=TRUE, col.names=TRUE,quote=F)
  write.table(confusion_matrix_and_indicators_arrange,paste(outfile_prefix, '_roc_confusion_matrix_and_indicators_arrange.txt', sep=""),sep="\t",row.names=TRUE, col.names=TRUE,quote=F)


  best_eval<-rownames(confusion_matrix_and_indicators_arrange)[str_detect(rownames(confusion_matrix_and_indicators_arrange),"_F1_score")][1]
  assign("best_eval_t_df",as.data.frame(t(as.data.frame(get(str_replace(best_eval,"_F1_score","_eval"))))))
  colnames(best_eval_t_df)<-"value_or_percentage"

  #show evaluation result as pdf table
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(fill = blues9[1:4], col=NA),
              fg_params=list(fontface=3)),
    colhead=list(fg_params=list(col="navyblue", fontface=4L)),
    rowhead=list(fg_params=list(col="orange", fontface=3L)))

  pdf(paste(outfile_prefix, '_roc_best_evaluation.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
  grid.arrange(
    tableGrob(best_eval_t_df, theme=tt3),
    nrow=1)
  invisible(dev.off())

  #read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
  evidence<-read.table(rRNAfile,head=F)#"hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed"
  colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

  #read hg38.psiU.SingleSites.bed
  evidence2<-read.table(rRNAfile2,head=F)#"hg38.psiU.SingleSites.bed"
  colnames(evidence2)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  evidence2$rRNA_uniq_id<-paste(evidence2$chrom,evidence2$chromStart,evidence2$chromEnd,evidence2$strand,sep="_")

  #filt by best F1 score threshold
  to_filt<-read.table(filtfile,head=F)
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
    "extendSeq")

  to_filt$roc_class<-as.numeric(to_filt[,str_replace(best_eval,"_F1_score","")]>get(str_replace(best_eval,"_F1_score","_thres")))
  to_filt$roc_class<-str_replace(as.character(to_filt$roc_class),"0","non-psi")
  to_filt$roc_class<-str_replace(as.character(to_filt$roc_class),"1","psi")
  table(to_filt$roc_class)
  write.table(to_filt,paste(outfile_prefix, '_roc_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  final_pred<-to_filt[to_filt[,str_replace(best_eval,"_F1_score","")]>get(str_replace(best_eval,"_F1_score","_thres")),]
  write.table(final_pred,paste(outfile_prefix, '_roc_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  write.table(final_pred,paste(outfile_prefix, '_roc_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

  #known_data miss/hit hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
  ROC_data$rRNA_uniq_id<-paste(ROC_data$chrom,ROC_data$chromStart,ROC_data$chromEnd,ROC_data$strand,sep="_")
  ROC_data_evidence<-ROC_data %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% as_tibble(.name_repair = "unique") %>% filter(!is.na(rRNA_anno))
  write.csv(ROC_data_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_hit.csv",sep=""))
  ROC_data_no_evidence<-evidence %>% left_join(ROC_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
  write.csv(ROC_data_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_miss.csv",sep=""))
  recovery<-paste(round(length(unique(ROC_data_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  #known_data miss/hit hg38.psiU.SingleSites.bed
  ROC_data$rRNA_uniq_id<-paste(ROC_data$chrom,ROC_data$chromStart,ROC_data$chromEnd,ROC_data$strand,sep="_")
  ROC_data_evidence2<-ROC_data %>% left_join(evidence2,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% as_tibble(.name_repair = "unique") %>% filter(!is.na(rRNA_anno))
  write.csv(ROC_data_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_hit.csv",sep=""))
  ROC_data_no_evidence2<-evidence2 %>% left_join(ROC_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
  write.csv(ROC_data_no_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_miss.csv",sep=""))
  recovery<-paste(round(length(unique(ROC_data_evidence2$rRNA_uniq_id))/length(unique(evidence2$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder recover (hg38.psiU.SingleSites.bed.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  #final_pred miss/hit
  final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
  final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% as_tibble(.name_repair = "unique") %>% filter(!is.na(rRNA_anno))
  write.csv(final_pred_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_roc_stopRatioFC_thres_hit.csv",sep=""))
  final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
  write.csv(final_pred_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_roc_stopRatioFC_thres_miss.csv",sep=""))
  recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiScan+ROC recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  print(paste("ROC evaluation in",output_dir,sep=" "))
}
