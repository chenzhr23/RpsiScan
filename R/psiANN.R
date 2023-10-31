# psiANN
#
# Build ANN model based on ground truth pseudouridylation data set

#' Build ANN model based on ground truth pseudouridylation data set
#' @param rocfile ROC input file of single sites information (with suffix _roc_plot.txt) generated from ground_truth() function
#' @param rRNAfile default we recommend hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
#' @param rRNAfile2 default we recommend hg38.psiU.SingleSites.bed
#' @param filtfile File to be filted
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
#'

#plot.nnet function from https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r
plot.nnet <- function(mod.in,nid=T,all.out=T,all.in=T,bias=T,wts.only=F,rel.rsc=5,circle.cex=5,
                      node.labs=T,var.labs=T,x.lab=NULL,y.lab=NULL,line.stag=NULL,struct=NULL,cex.val=1,
                      alpha.val=1,circle.col='lightblue',pos.col='black',neg.col='grey', max.sp = F, ...){

  #sanity checks
  if('mlp' %in% class(mod.in)) warning('Bias layer not applicable for rsnns object')
  if('numeric' %in% class(mod.in)){
    if(is.null(struct)) stop('Three-element vector required for struct')
    if(length(mod.in) != ((struct[1]*struct[2]+struct[2]*struct[3])+(struct[3]+struct[2])))
      stop('Incorrect length of weight matrix for given network structure')
  }
  if('train' %in% class(mod.in)){
    if('nnet' %in% class(mod.in$finalModel)){
      mod.in<-mod.in$finalModel
      warning('Using best nnet model from train output')
    }
    else stop('Only nnet method can be used with train object')
  }

  #gets weights for neural network, output is list
  #if rescaled argument is true, weights are returned but rescaled based on abs value
  nnet.vals<-function(mod.in,nid,rel.rsc,struct.out=struct){

    if('numeric' %in% class(mod.in)){
      struct.out<-struct
      wts<-mod.in
    }

    #neuralnet package
    if('nn' %in% class(mod.in)){
      struct.out<-unlist(lapply(mod.in$weights[[1]],ncol))
      struct.out<-struct.out[-length(struct.out)]
      struct.out<-c(
        length(mod.in$model.list$variables),
        struct.out,
        length(mod.in$model.list$response)
      )
      wts<-unlist(mod.in$weights[[1]])
    }

    #nnet package
    if('nnet' %in% class(mod.in)){
      struct.out<-mod.in$n
      wts<-mod.in$wts
    }

    #RSNNS package
    if('mlp' %in% class(mod.in)){
      struct.out<-c(mod.in$nInputs,mod.in$archParams$size,mod.in$nOutputs)
      hid.num<-length(struct.out)-2
      wts<-mod.in$snnsObject$getCompleteWeightMatrix()

      #get all input-hidden and hidden-hidden wts
      inps<-wts[grep('Input',row.names(wts)),grep('Hidden_2',colnames(wts)),drop=F]
      inps<-melt(rbind(rep(NA,ncol(inps)),inps))$value
      uni.hids<-paste0('Hidden_',1+seq(1,hid.num))
      for(i in 1:length(uni.hids)){
        if(is.na(uni.hids[i+1])) break
        tmp<-wts[grep(uni.hids[i],rownames(wts)),grep(uni.hids[i+1],colnames(wts)),drop=F]
        inps<-c(inps,melt(rbind(rep(NA,ncol(tmp)),tmp))$value)
      }

      #get connections from last hidden to output layers
      outs<-wts[grep(paste0('Hidden_',hid.num+1),row.names(wts)),grep('Output',colnames(wts)),drop=F]
      outs<-rbind(rep(NA,ncol(outs)),outs)

      #weight vector for all
      wts<-c(inps,melt(outs)$value)
      assign('bias',F,envir=environment(nnet.vals))
    }

    if(nid) wts<-rescale(abs(wts),c(1,rel.rsc))

    #convert wts to list with appropriate names
    hid.struct<-struct.out[-c(length(struct.out))]
    row.nms<-NULL
    for(i in 1:length(hid.struct)){
      if(is.na(hid.struct[i+1])) break
      row.nms<-c(row.nms,rep(paste('hidden',i,seq(1:hid.struct[i+1])),each=1+hid.struct[i]))
    }
    row.nms<-c(
      row.nms,
      rep(paste('out',seq(1:struct.out[length(struct.out)])),each=1+struct.out[length(struct.out)-1])
    )
    out.ls<-data.frame(wts,row.nms)
    out.ls$row.nms<-factor(row.nms,levels=unique(row.nms),labels=unique(row.nms))
    out.ls<-split(out.ls$wts,f=out.ls$row.nms)

    assign('struct',struct.out,envir=environment(nnet.vals))

    out.ls

  }

  wts<-nnet.vals(mod.in,nid=F)

  if(wts.only) return(wts)

  #circle colors for input, if desired, must be two-vector list, first vector is for input layer
  if(is.list(circle.col)){
    circle.col.inp<-circle.col[[1]]
    circle.col<-circle.col[[2]]
  }
  else circle.col.inp<-circle.col

  #initiate plotting
  x.range<-c(0,100)
  y.range<-c(0,100)
  #these are all proportions from 0-1
  if(is.null(line.stag)) line.stag<-0.011*circle.cex/2
  layer.x<-seq(0.17,0.9,length=length(struct))
  bias.x<-layer.x[-length(layer.x)]+diff(layer.x)/2
  bias.y<-0.95
  circle.cex<-circle.cex

  #get variable names from mod.in object
  #change to user input if supplied
  if('numeric' %in% class(mod.in)){
    x.names<-paste0(rep('X',struct[1]),seq(1:struct[1]))
    y.names<-paste0(rep('Y',struct[3]),seq(1:struct[3]))
  }
  if('mlp' %in% class(mod.in)){
    all.names<-mod.in$snnsObject$getUnitDefinitions()
    x.names<-all.names[grep('Input',all.names$unitName),'unitName']
    y.names<-all.names[grep('Output',all.names$unitName),'unitName']
  }
  if('nn' %in% class(mod.in)){
    x.names<-mod.in$model.list$variables
    y.names<-mod.in$model.list$respons
  }
  if('xNames' %in% names(mod.in)){
    x.names<-mod.in$xNames
    y.names<-attr(terms(mod.in),'factor')
    y.names<-row.names(y.names)[!row.names(y.names) %in% x.names]
  }
  if(!'xNames' %in% names(mod.in) & 'nnet' %in% class(mod.in)){
    if(is.null(mod.in$call$formula)){
      x.names<-colnames(eval(mod.in$call$x))
      y.names<-colnames(eval(mod.in$call$y))
    }
    else{
      forms<-eval(mod.in$call$formula)
      x.names<-mod.in$coefnames
      facts<-attr(terms(mod.in),'factors')
      y.check<-mod.in$fitted
      if(ncol(y.check)>1) y.names<-colnames(y.check)
      else y.names<-as.character(forms)[2]
    }
  }
  #change variables names to user sub
  if(!is.null(x.lab)){
    if(length(x.names) != length(x.lab)) stop('x.lab length not equal to number of input variables')
    else x.names<-x.lab
  }
  if(!is.null(y.lab)){
    if(length(y.names) != length(y.lab)) stop('y.lab length not equal to number of output variables')
    else y.names<-y.lab
  }

  #initiate plot
  plot(x.range,y.range,type='n',axes=F,ylab='',xlab='',...)

  #function for getting y locations for input, hidden, output layers
  #input is integer value from 'struct'
  get.ys<-function(lyr, max_space = max.sp){
    if(max_space){
      spacing <- diff(c(0*diff(y.range),0.9*diff(y.range)))/lyr
    } else {
      spacing<-diff(c(0*diff(y.range),0.9*diff(y.range)))/max(struct)
    }

    seq(0.5*(diff(y.range)+spacing*(lyr-1)),0.5*(diff(y.range)-spacing*(lyr-1)),
        length=lyr)
  }

  #function for plotting nodes
  #'layer' specifies which layer, integer from 'struct'
  #'x.loc' indicates x location for layer, integer from 'layer.x'
  #'layer.name' is string indicating text to put in node
  layer.points<-function(layer,x.loc,layer.name,cex=cex.val){
    x<-rep(x.loc*diff(x.range),layer)
    y<-get.ys(layer)
    points(x,y,pch=21,cex=circle.cex,col=in.col,bg=bord.col)
    if(node.labs) text(x,y,paste(layer.name,1:layer,sep=''),cex=cex.val)
    if(layer.name=='I' & var.labs) text(x-line.stag*diff(x.range),y,x.names,pos=2,cex=cex.val)
    if(layer.name=='O' & var.labs) text(x+line.stag*diff(x.range),y,y.names,pos=4,cex=cex.val)
  }

  #function for plotting bias points
  #'bias.x' is vector of values for x locations
  #'bias.y' is vector for y location
  #'layer.name' is  string indicating text to put in node
  bias.points<-function(bias.x,bias.y,layer.name,cex,...){
    for(val in 1:length(bias.x)){
      points(
        diff(x.range)*bias.x[val],
        bias.y*diff(y.range),
        pch=21,col=in.col,bg=bord.col,cex=circle.cex
      )
      if(node.labs)
        text(
          diff(x.range)*bias.x[val],
          bias.y*diff(y.range),
          paste(layer.name,val,sep=''),
          cex=cex.val
        )
    }
  }

  #function creates lines colored by direction and width as proportion of magnitude
  #use 'all.in' argument if you want to plot connection lines for only a single input node
  layer.lines<-function(mod.in,h.layer,layer1=1,layer2=2,out.layer=F,nid,rel.rsc,all.in,pos.col,
                        neg.col,...){

    x0<-rep(layer.x[layer1]*diff(x.range)+line.stag*diff(x.range),struct[layer1])
    x1<-rep(layer.x[layer2]*diff(x.range)-line.stag*diff(x.range),struct[layer1])

    if(out.layer==T){

      y0<-get.ys(struct[layer1])
      y1<-rep(get.ys(struct[layer2])[h.layer],struct[layer1])
      src.str<-paste('out',h.layer)

      wts<-nnet.vals(mod.in,nid=F,rel.rsc)
      wts<-wts[grep(src.str,names(wts))][[1]][-1]
      wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
      wts.rs<-wts.rs[grep(src.str,names(wts.rs))][[1]][-1]

      cols<-rep(pos.col,struct[layer1])
      cols[wts<0]<-neg.col

      if(nid) segments(x0,y0,x1,y1,col=cols,lwd=wts.rs)
      else segments(x0,y0,x1,y1)

    }

    else{

      if(is.logical(all.in)) all.in<-h.layer
      else all.in<-which(x.names==all.in)

      y0<-rep(get.ys(struct[layer1])[all.in],struct[2])
      y1<-get.ys(struct[layer2])
      src.str<-paste('hidden',layer1)

      wts<-nnet.vals(mod.in,nid=F,rel.rsc)
      wts<-unlist(lapply(wts[grep(src.str,names(wts))],function(x) x[all.in+1]))
      wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)
      wts.rs<-unlist(lapply(wts.rs[grep(src.str,names(wts.rs))],function(x) x[all.in+1]))

      cols<-rep(pos.col,struct[layer2])
      cols[wts<0]<-neg.col

      if(nid) segments(x0,y0,x1,y1,col=cols,lwd=wts.rs)
      else segments(x0,y0,x1,y1)

    }

  }

  bias.lines<-function(bias.x,mod.in,nid,rel.rsc,all.out,pos.col,neg.col,...){

    if(is.logical(all.out)) all.out<-1:struct[length(struct)]
    else all.out<-which(y.names==all.out)

    for(val in 1:length(bias.x)){

      wts<-nnet.vals(mod.in,nid=F,rel.rsc)
      wts.rs<-nnet.vals(mod.in,nid=T,rel.rsc)

      if(val != length(bias.x)){
        wts<-wts[grep('out',names(wts),invert=T)]
        wts.rs<-wts.rs[grep('out',names(wts.rs),invert=T)]
        sel.val<-grep(val,substr(names(wts.rs),8,8))
        wts<-wts[sel.val]
        wts.rs<-wts.rs[sel.val]
      }

      else{
        wts<-wts[grep('out',names(wts))]
        wts.rs<-wts.rs[grep('out',names(wts.rs))]
      }

      cols<-rep(pos.col,length(wts))
      cols[unlist(lapply(wts,function(x) x[1]))<0]<-neg.col
      wts.rs<-unlist(lapply(wts.rs,function(x) x[1]))

      if(nid==F){
        wts.rs<-rep(1,struct[val+1])
        cols<-rep('black',struct[val+1])
      }

      if(val != length(bias.x)){
        segments(
          rep(diff(x.range)*bias.x[val]+diff(x.range)*line.stag,struct[val+1]),
          rep(bias.y*diff(y.range),struct[val+1]),
          rep(diff(x.range)*layer.x[val+1]-diff(x.range)*line.stag,struct[val+1]),
          get.ys(struct[val+1]),
          lwd=wts.rs,
          col=cols
        )
      }

      else{
        segments(
          rep(diff(x.range)*bias.x[val]+diff(x.range)*line.stag,struct[val+1]),
          rep(bias.y*diff(y.range),struct[val+1]),
          rep(diff(x.range)*layer.x[val+1]-diff(x.range)*line.stag,struct[val+1]),
          get.ys(struct[val+1])[all.out],
          lwd=wts.rs[all.out],
          col=cols[all.out]
        )
      }

    }
  }

  #use functions to plot connections between layers
  #bias lines
  if(bias) bias.lines(bias.x,mod.in,nid=nid,rel.rsc=rel.rsc,all.out=all.out,pos.col=alpha(pos.col,alpha.val),
                      neg.col=alpha(neg.col,alpha.val))

  #layer lines, makes use of arguments to plot all or for individual layers
  #starts with input-hidden
  #uses 'all.in' argument to plot connection lines for all input nodes or a single node
  if(is.logical(all.in)){
    mapply(
      function(x) layer.lines(mod.in,x,layer1=1,layer2=2,nid=nid,rel.rsc=rel.rsc,
                              all.in=all.in,pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val)),
      1:struct[1]
    )
  }
  else{
    node.in<-which(x.names==all.in)
    layer.lines(mod.in,node.in,layer1=1,layer2=2,nid=nid,rel.rsc=rel.rsc,all.in=all.in,
                pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val))
  }
  #connections between hidden layers
  lays<-split(c(1,rep(2:(length(struct)-1),each=2),length(struct)),
              f=rep(1:(length(struct)-1),each=2))
  lays<-lays[-c(1,(length(struct)-1))]
  for(lay in lays){
    for(node in 1:struct[lay[1]]){
      layer.lines(mod.in,node,layer1=lay[1],layer2=lay[2],nid=nid,rel.rsc=rel.rsc,all.in=T,
                  pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val))
    }
  }
  #lines for hidden-output
  #uses 'all.out' argument to plot connection lines for all output nodes or a single node
  if(is.logical(all.out))
    mapply(
      function(x) layer.lines(mod.in,x,layer1=length(struct)-1,layer2=length(struct),out.layer=T,nid=nid,rel.rsc=rel.rsc,
                              all.in=all.in,pos.col=alpha(pos.col,alpha.val),neg.col=alpha(neg.col,alpha.val)),
      1:struct[length(struct)]
    )
  else{
    node.in<-which(y.names==all.out)
    layer.lines(mod.in,node.in,layer1=length(struct)-1,layer2=length(struct),out.layer=T,nid=nid,rel.rsc=rel.rsc,
                pos.col=pos.col,neg.col=neg.col,all.out=all.out)
  }

  #use functions to plot nodes
  for(i in 1:length(struct)){
    in.col<-bord.col<-circle.col
    layer.name<-'H'
    if(i==1) { layer.name<-'I'; in.col<-bord.col<-circle.col.inp}
    if(i==length(struct)) layer.name<-'O'
    layer.points(struct[i],layer.x[i],layer.name)
  }

  if(bias) bias.points(bias.x,bias.y,'B')

}

psiANN <- function(rocfile, rRNAfile, rRNAfile2, filtfile,output_dir,output_name)
{
  #nohup Rscript $(dirname "$0")/ann_build_totalRNA.r -f ${outFileName}_ann_plot.txt -k ${outFileName}_ann_filt_totalRNA.bed -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed  -o ${outFileName} > ${outFileName}_ann_evaluation_totalRNA.log 2>&1 &

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  outfile_prefix<-paste(output_dir,"/",output_name,sep="")

  #load data
  cat("\n\n=====================Load data to perform classification model (factor as input)=====================\n")
  ANN_data<-as.data.frame(fread(rocfile))
  colnames(ANN_data)<-c(
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
  # ANN_data$base<-"T"
  ANN_data$base<-str_replace_all(ANN_data$base, "TRUE", "T")

  ANN_data_sel<-ANN_data %>% select(tretRtsRatio,rtsRatioFold,tretBefRatio,befRatioFold,tretAftRatio,aftRatioFold,tretMutRatio,mutRatioFold,tretDelRatio,delRatioFold,class)
  rownames(ANN_data_sel)<-ANN_data$name

  ANN_data_sel$class<-str_replace_all(ANN_data_sel$class, "0", "non_psi")
  ANN_data_sel$class<-str_replace_all(ANN_data_sel$class, "1", "psi")

  #Encoding the target feature as factor
  ANN_data_sel$class<-as.factor(ANN_data_sel$class)
  cat("total classification: ")
  table(ANN_data_sel$class)


  # Splitting the dataset into the Training set and Test set
  set.seed(123)
  ANN_data_sel$psi <- ANN_data_sel$class == "psi"
  ANN_data_sel$non_psi <- ANN_data_sel$class == "non_psi"
  split = sample.split(ANN_data_sel$class, SplitRatio = 0.7)
  training_set = subset(ANN_data_sel, split == TRUE)
  test_set = subset(ANN_data_sel, split == FALSE)

  # feature scaling
  training_set_mean <- apply(training_set %>% select(-class, -psi, -non_psi), 2, mean)
  training_set_sd <- apply(training_set %>% select(-class, -psi, -non_psi), 2, sd)
  training_set[, c(-11, -12, -13)] = scale(training_set[, c(-11, -12, -13)], center = training_set_mean, scale = training_set_sd)
  test_set[, c(-11, -12, -13)] = scale(test_set[, c(-11, -12, -13)], center = training_set_mean, scale = training_set_sd)

  training_set_origin = subset(ANN_data, split == TRUE)
  test_set_origin = subset(ANN_data, split == FALSE)
  cat("\n", "training set classification using sample.split: ")
  table(training_set$class)
  summary(training_set)
  cat("\n", "test set classification using sample.split: ")
  table(test_set$class)
  summary(test_set)

  #Network Aplication
  # ANN_data_sel.net <- neuralnet(psi + non_psi ~
  #                       treatPreRpmFold + preRpmFoldRatio + treatAfterRpmFold + afterRpmFoldRatio + treatStopRatio + stopRatioFC,
  #                       data = training_set, hidden = c(12, 12), rep = 5, err.fct = "ce",
  #                       linear.output = FALSE, lifesign = "minimal", stepmax = 1000000,
  #                       threshold = 0.001)

  ANN_data_sel.net <- neuralnet(psi + non_psi ~ tretRtsRatio + rtsRatioFold + tretBefRatio + befRatioFold + tretAftRatio + aftRatioFold + tretMutRatio + mutRatioFold + tretDelRatio + delRatioFold,data = training_set, hidden = c(12, 12, 12, 12), rep = 10)

  #Get weights for a neural network in an organized list by extracting values from a neural network object
  neuralweights(ANN_data_sel.net)

  #import the function from Github
  # source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
  #plot each model
  pdf(paste(outfile_prefix, '_ann_nnet_train_data_plot.pdf', sep = ""),width=16,height=8)
  plot.nnet(ANN_data_sel.net)
  invisible(dev.off())

  #Visualize the network structure of the fully connected model
  pdf(paste(outfile_prefix, '_ann_best_train_data_plot.pdf', sep = ""))
  plot(ANN_data_sel.net, rep = "best")
  invisible(dev.off())

  #Visualize the prediction effect of the model
  mlppre <- predict(ANN_data_sel.net, test_set)
  colnames(mlppre)<-c("psi","non_psi")
  mlpprelab <- apply(mlppre, 1, which.max)
  mlpprelab<-sub("1","psi",as.character(mlpprelab))
  mlpprelab<-sub("2","non_psi",as.character(mlpprelab))
  plot_confusion_matrix <- ggplot() +
    geom_confmat(aes(x = test_set$class, y = mlpprelab),
                 normalize = TRUE, text.perc = TRUE) +
    labs(x = "Reference", y = "Prediction") +
    scale_fill_gradient2(low = "darkblue", high = "#ec2f2f") +
    theme_bw() +
    theme(plot.margin = unit(c(6, 5, 6, 5), "cm"))

  pdf(paste(outfile_prefix, '_ann_best_test_confusion_matrix.pdf', sep = ""))
  plot_confusion_matrix
  invisible(dev.off())

  # Predicting Result
  ANN_data_sel.prediction <- neuralnet::compute(ANN_data_sel.net, test_set[-11:-13])
  idx <- apply(ANN_data_sel.prediction$net.result, 1, which.max)
  net.result <- as.data.frame(ANN_data_sel.prediction$net.result)
  colnames(net.result) <- c('psi', 'non_psi')
  pred <- c('psi', 'non_psi')[idx]
  table(pred)

  pdf(paste(outfile_prefix, '_ann_roc_test_data_plot.pdf', sep=""))
  ann_roc_test <- roc(main = "Test Data ROC", test_set$class, net.result$psi, smooth = FALSE, print.auc = TRUE, col = "#e41a1c", plot = TRUE, print.thres = "best", print.thres.best.method = "youden", levels = c("non_psi", "psi"), direction = '<', auc = T, ci = T)
  invisible(dev.off())
  cat("best roc threshold: ")
  coords(ann_roc_test, "best", ret = "all", transpose = TRUE)[1]

  summary(ANN_data_sel.net)
  str(ANN_data_sel.net)
  save(ANN_data_sel.net, ann_roc_test, training_set_mean, training_set_sd, file = paste(outfile_prefix, '_ANN_model.RData', sep = "")) #"my-nn.RData"
  # saveRDS(mymodel, file = paste(outfile_prefix, '_ann_model.rds', sep=""))#"my-nn.rds"

  #add attr to ann_test_data
  ann_test_data <- test_set_origin
  # ann_test_data$ann_test_prob_class <- ifelse(net.result$psi >= coords(ann_roc_test, "best", ret = "all", transpose = TRUE)[1], "psi", "non_psi")
  ann_test_data <- cbind(ann_test_data, net.result, pred)
  write.xlsx(ann_test_data, paste(outfile_prefix, '_ann_test_data.xlsx', sep = ""), overwrite = TRUE)

  #output model info
  cat("\n\n=====================tune best model confusion matrix (for modeling)=====================\n")
  tab <- table(Predicted = pred, Actual = test_set$class)
  tab
  cat("tune best model error rate (for modeling): ", 1 - sum(diag(tab)) / sum(tab), "\n")
  cat("tune best model correct rate (for modeling): ", sum(diag(tab)) / sum(tab), "\n")

  #calculate evaluation indicators
  cat("\n\n=====================Calculate evaluation indicators=====================\n")
  confusion_matrix <- as.data.frame(tab)
  ann_TP <- confusion_matrix[confusion_matrix$Predicted == "psi" & confusion_matrix$Actual == "psi",]$Freq #True Positives (TP)
  ann_FP <- confusion_matrix[confusion_matrix$Predicted == "psi" & confusion_matrix$Actual == "non_psi",]$Freq #False Positives (FP)
  ann_TN <- confusion_matrix[confusion_matrix$Predicted == "non_psi" & confusion_matrix$Actual == "non_psi",]$Freq #True Negatives (TN)
  ann_FN <- confusion_matrix[confusion_matrix$Predicted == "non_psi" & confusion_matrix$Actual == "psi",]$Freq #False Negatives (FN)
  ann_TPR <- ann_TP / (ann_TP + ann_FN) #sensitivity (true positive rate, TPR)
  ann_TNR <- ann_TN / (ann_TN + ann_FP) #specifity (selectivity or true negative rate, TNR)
  ann_FPR <- 1 - ann_TNR #False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
  ann_FNR <- 1 - ann_TPR #False Negative Rate, FNR)
  ann_Prec <- ann_TP / (ann_TP + ann_FP) #Precision
  ann_Recall <- ann_TP / (ann_TP + ann_FN) #Recall
  ann_ACC <- (ann_TP + ann_TN) / (ann_TP + ann_TN + ann_FP + ann_FN) #accuracy
  ann_F1_score <- (2 * ann_Recall * ann_Prec) / (ann_Recall + ann_Prec) #F1_score
  eval <- cbind(ann_TP, ann_FP, ann_TN, ann_FN, ann_TPR, ann_TNR, ann_FPR, ann_FNR, ann_Prec, ann_Recall, ann_ACC, ann_F1_score)
  eval <- round(eval, 3)
  eval
  write.table(eval, paste(outfile_prefix, '_ann_eval.txt', sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)

  #show nn evaluation as pdf table
  ann_eval_t_df <- as.data.frame(t(as.data.frame(eval)))
  colnames(ann_eval_t_df) <- "value_or_percentage"
  tt3 <- ttheme_minimal(
    core = list(bg_params = list(fill = blues9[1:4], col = NA),
                fg_params = list(fontface = 3)),
    colhead = list(fg_params = list(col = "navyblue", fontface = 4L)),
    rowhead = list(fg_params = list(col = "orange", fontface = 3L)))

  pdf(paste(outfile_prefix, '_ann_evaluation.pdf', sep = ""), width = 7, height = 7) # Open a new pdf file
  grid.arrange(
    tableGrob(ann_eval_t_df, theme = tt3),
    nrow = 1)
  invisible(dev.off()) # Close the file


  #filt by nn best model
  cat("\n\n=====================Filt by nn best model=====================\n")
  cat("below is your input data ready to be predicted...\n")
  # get prediction
  to_pred <- as.data.frame(fread(filtfile))
  colnames(to_pred)<-c(
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
  to_pred$base<-str_replace_all(to_pred$base, "TRUE", "T")

  to_pred_var <- to_pred %>% select(tretRtsRatio,rtsRatioFold,tretBefRatio,befRatioFold,tretAftRatio,aftRatioFold,tretMutRatio,mutRatioFold,tretDelRatio,delRatioFold)
  to_pred_var = scale(to_pred_var, center = training_set_mean, scale = training_set_sd)
  str(to_pred_var)

  # Predicting Result
  ANN_data_sel.prediction <- neuralnet::compute(ANN_data_sel.net, to_pred_var)
  idx <- apply(ANN_data_sel.prediction$net.result, 1, which.max)
  net.result <- as.data.frame(ANN_data_sel.prediction$net.result)
  colnames(net.result) <- c('psi', 'non_psi')
  pred <- c('psi', 'non_psi')[idx]
  table(pred)


  #add attr to ann_pred_data
  ann_pred_data <- to_pred
  ann_pred_data <- cbind(ann_pred_data, net.result, pred)
  write.xlsx(ann_pred_data, paste(outfile_prefix, '_ann_total_prediction.xlsx', sep = ""), overwrite = TRUE)

  #read training/test set evidence
  evidence<-read.delim(rRNAfile,head=F)#"training/test set"
  colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

  #read hg38.psiU.SingleSites.bed
  evidence2<-read.delim(rRNAfile2,head=F)#"hg38.psiU.SingleSites.bed"
  colnames(evidence2)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  evidence2$rRNA_uniq_id<-paste(evidence2$chrom,evidence2$chromStart,evidence2$chromEnd,evidence2$strand,sep="_")

  #known_data miss/hit hg38_human_chr21_rRNA_known_pseudoU_SingleSites
  ANN_data$rRNA_uniq_id<-paste(ANN_data$chrom,ANN_data$chromStart,ANN_data$chromEnd,ANN_data$strand,sep="_")
  ANN_data_evidence<-ANN_data %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
  write.csv(ANN_data_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_hit.csv",sep=""))
  ANN_data_no_evidence<-evidence %>% left_join(ANN_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
  write.csv(ANN_data_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_miss.csv",sep=""))
  recovery<-paste(round(length(unique(ANN_data_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder recover (training/test set)",recovery,"rRNA psi sites in all known chrom21\n")

  #known_data miss/hit hg38.psiU.SingleSites.bed
  ANN_data$rRNA_uniq_id<-paste(ANN_data$chrom,ANN_data$chromStart,ANN_data$chromEnd,ANN_data$strand,sep="_")
  ANN_data_evidence2<-ANN_data %>% left_join(evidence2,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
  write.csv(ANN_data_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_hit.csv",sep=""))
  ANN_data_no_evidence2<-evidence2 %>% left_join(ANN_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
  write.csv(ANN_data_no_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_miss.csv",sep=""))
  recovery<-paste(round(length(unique(ANN_data_evidence2$rRNA_uniq_id))/length(unique(evidence2$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder recover (hg38.psiU.SingleSites.bed.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  #final_pred miss/hit
  final_pred<-ann_pred_data[ann_pred_data$pred=="psi",]
  final_pred<-final_pred[final_pred$foldChange>2,]
  write.table(final_pred,paste(outfile_prefix, '_ann_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  write.table(final_pred,paste(outfile_prefix, '_ann_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

  final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
  final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
  write.csv(final_pred_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_ann_hit.csv",sep=""))
  final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
  write.csv(final_pred_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_ann_miss.csv",sep=""))
  recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder+ANN recover (training/test set)",recovery,"rRNA psi sites in all known chrom21\n")

  print(paste("ANN model in",output_dir,sep=" "))
}
