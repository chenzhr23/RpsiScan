# RpsiScan (CIV-seq)

## Solution to Ψ-site identification, annotation, and target prediction

## Contents
- [pre-installation](#pre-installation)
- [test data](#test-data)
- [Installation](#Installation)
- [Usage](#Usage)
  - [Ψ-site identification](#Ψ-site-identification)
    - [calculate overall reverse transcription stop information](#calculate-overall-reverse-transcription-stop-information)
    - [pre-filt overall reverse transcription stop information](#pre-filt-overall-reverse-transcription-stop-information)
    - [a) User-defined](#a-User-defined)
    - [b) SVM](#b-SVM)
    - [c) ANN](#c-ANN)
  - [Ψ-site annotation](#Ψ-site-annotation)
  - [Ψ-site target prediction](#Ψ-site-target-prediction)

### pre-installation
RpsiScan R package is a R API mostly for psiScan, psiAnnotator, and psiACAscan (also include bedtools utilization for some functions). psiScan/psiAnnotator/psiACAscan are C programs and predominantly used in unix-based operating systems. Therefore, for the usability of RpsiScan, we recommend install RpsiScan R package in WSL2 (WSL2 installation guide: https://pureinfotech.com/install-windows-subsystem-linux-2-windows-10/).

### test data
Test data can be downloaded from https://mega.nz/fm/public-links/YGVghDRA.

### Installation

```R
#use pacman to install packages in batch

install.packages("pacman")
library(pacman)

#load and install required R packages 
p_load('remotes','dplyr','data.table','pROC','ggplot2','cowplot','gridExtra','reshape2','stringr','RColorBrewer','scales','e1071','openxlsx','caTools','ggpol','neuralnet','NeuralNetTools','tidyr')

#install from github
library("remotes")#or run library("devtools")
install_github("chenzhr23/RpsiScan")
```

### Usage

#### Ψ-site identification

##### calculate overall reverse transcription stop information

RpsiScan function was used to generate overall RT site information with the following default options: -p 1.5 -t 5 -r 0.05 -M 1 -f 1 -m 0 -s -w 20. gene_bed12: bed12 file of annotation, usually downloaded form databases like UCSC, ENSEMBL, etc., or organize the file by bioinformatics tools

```R
#run RpsiScan and invoke psiScan C program to calculate overall reverse transcription stop information
#RpsiScan will generate a file in output_dir/ and named as output_name.txt (for example: ../test_out/RpsiScan_out/RpsiScan_out.txt)
#when RpsiScan function is done, the global environment will assign a variable called 'output_name' (for example: RpsiScan_out), which contains the overall reverse transcription stop information (RTS information derived from CMC-treated/CMC-input dataset)
library("RpsiScan")
RpsiScan(genome_fa="../test_data/genome/hg38.fa",
            genome_fai="../test_data/genome/hg38.fa.fai",
            treat_bam="../test_data/C2.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gzAligned.sortedByCoord.out.bam",
            input_bam = "../test_data/C1.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gzAligned.sortedByCoord.out.bam",
            output_dir="../test_out/RpsiScan_out",
            output_name = "RpsiScan_out")

#or load the demoe data we deposit in the package
RpsiScan_example()#when finished, view your global environment, and you will detect Data named 'RpsiScan_example' (this will load ePSI-seq_test.txt, RTS information from a ePSI-seq total RNA sample)
```

##### pre-filt overall reverse transcription stop information

Remove sites with low RT stop read abundance to avoid false positives caused by too small RT reads number (treatStopNum > 10); only retains called sites of U base to include the enrichment sites caused by Ψ-CMC (not by G-CMC or others)

```R
#prefilt function will retain psiScan information with desired target chrom/target base.
prefilt_out <- prefilt(RpsiScan_res=RpsiScan_example,
                        target_chrom="^chr[0-9|a-z|A-Z]*$",
                        target_base="T",
                        output_dir="../test_out/RpsiScan_out",
                        output_name="test_prefilt")
```
Once get the prefilt result, users can choose one of the approach to perform Ψ-site identification: a) User-defined, b) SVM, c) ANN.

##### a) User-defined

```R
RpsiScan_example()#when finished, view your global environment, and you will detect Data named 'RpsiScan_example'

prefilt_out <- prefilt(RpsiScan_res=RpsiScan_example,
                        target_chrom="^chr[0-9|a-z|A-Z]*$",
                        target_base="T",
                        output_dir="../test_out/RpsiScan_out",
                        output_name="test_prefilt")

#udfilt function will retain psiScan information with desired identification metrics greater than thresholds (treatPreRpmFold, treatAfterRpmFold, treatStopRatio, preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC).
udfilt_out <- udfilt(RpsiScan_res_file="../test_out/RpsiScan_out/test_prefilt.txt",
                   tretRtsRatio_thres=0.02,
                   rtsRatioFold_thres=1.5,
                   tretBefRatio_thres=0.02,
                   befRatioFold_thres=1.5,
                   tretAftRatio_thres=0.02,
                   aftRatioFold_thres=1.5,
                   tretMutRatio_thres=0.02,
                   mutRatioFold_thres=1.5,
                   tretDelRatio_thres=0.02,
                   delRatioFold_thres=1.5,
                   output_dir="../test_out/RpsiScan_out",
                   output_name="test_udfilt")

#ground_truth function prepare ground truth data set (known rRNA positive/negative pseudouridylation sites) for pseudouridylation identification.
ground_truth(RpsiScan_res_file="../test_out/RpsiScan_out/test_udfilt.txt",
              known_rRNA_bed="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
              rrna_chr21_bed="../test_data/ground_truth/rrna_chr21.bed",
              output_dir="../test_out/ground_truth",
              output_name="ground_truth")

#user_defined function will accept result file generated by ground_truth() function and give out evaluation on input ground truth rRNA pseudouridylation data set.
user_defined(rocfile="../test_out/ground_truth/ground_truth_roc_plot.txt",
              rRNAfile="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
              filtfile="../test_out/RpsiScan_out/test_udfilt.txt",
              output_dir="../test_out/user_defined",
              output_name="user_defined")
```

##### b) SVM

```R
RpsiScan_example()#when finished, view your global environment, and you will detect Data named 'RpsiScan_example'

prefilt_out <- prefilt(RpsiScan_res=RpsiScan_example,
                        target_chrom="^chr[0-9|a-z|A-Z]*$",
                        target_base="T",
                        output_dir="../test_out/RpsiScan_out",
                        output_name="test_prefilt")

ground_truth(RpsiScan_res_file="../test_out/RpsiScan_out/test_prefilt.txt",
              known_rRNA_bed="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
              rrna_chr21_bed="../test_data/ground_truth/rrna_chr21.bed",
              output_dir="../test_out/ground_truth",
              output_name="ground_truth")

#Perform ROC evaluation based on ground truth pseudouridylation data set
psiROC(rocfile="../test_out/ground_truth/ground_truth_roc_plot.txt",
        rRNAfile="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
        rRNAfile2="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
        filtfile="../test_out/RpsiScan_out/test_prefilt.txt",
        output_dir="../test_out/roc",
        output_name="roc")

#Build SVM model based on ground truth pseudouridylation data set
psiSVM(rocfile="../test_out/ground_truth/ground_truth_roc_plot.txt",
        rRNAfile="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
        rRNAfile2="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
        filtfile="../test_out/RpsiScan_out/test_prefilt.txt",
        output_dir="../test_out/svm_model",
        output_name="svm_model")

#Predict pseudouridylation sites based on SVM model
psiSVMpred(filtfile="../test_out/RpsiScan_out/test_prefilt.txt",
           svmmodelfile="../test_out/svm_model/svm_model_SVM_model.RData",
           rRNAfile="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
           output_dir="../test_out/svm_predict",
           output_name="svm_predict")

```

##### c) ANN

```R
RpsiScan_example()#when finished, view your global environment, and you will detect Data named 'RpsiScan_example'

prefilt_out <- prefilt(RpsiScan_res=RpsiScan_example,
                        target_chrom="^chr[0-9|a-z|A-Z]*$",
                        target_base="T",
                        output_dir="../test_out/RpsiScan_out",
                        output_name="test_prefilt")

ground_truth(RpsiScan_res_file="../test_out/RpsiScan_out/test_prefilt.txt",
              known_rRNA_bed="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
              rrna_chr21_bed="../test_data/ground_truth/rrna_chr21.bed",
              output_dir="../test_out/ground_truth",
              output_name="ground_truth")

#Perform ROC evaluation based on ground truth pseudouridylation data set
psiROC(rocfile="../test_out/ground_truth/ground_truth_roc_plot.txt",
        rRNAfile="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
        rRNAfile2="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
        filtfile="../test_out/RpsiScan_out/test_prefilt.txt",
        output_dir="../test_out/roc",
        output_name="roc")

#Build ANN model based on ground truth pseudouridylation data set
psiANN(rocfile="../test_out/ground_truth/ground_truth_roc_plot.txt",
        rRNAfile="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
        rRNAfile2="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
        filtfile="../test_out/RpsiScan_out/test_prefilt.txt",
        output_dir="../test_out/ann_model",
        output_name="ann_model")

#Predict pseudouridylation sites based on ANN model
psiANNpred(filtfile="../test_out/RpsiScan_out/test_prefilt.txt",
           annmodelfile="../test_out/ann_model/ann_model_ANN_model.RData",
           rRNAfile="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
           output_dir="../test_out/ann_predict",
           output_name="ann_predict")
```

#### Ψ-site annotation
```R
#annotate psiScan result generated by User-defined thresholds
RpsiAnnotator(input_bedfile="../test_out/RpsiScan_out/test_udfilt.bed",
            annotation_bed6file="../test_data/annotation/hg38.genecode.v30.tRNA.snoRNA.miRNA.rmsk.exonFeatures.bed6",
            output_dir="../test_out/ud_annotation",
            output_name = "ud_annotation")

#annotate prediction result generated by ROC best threshold
RpsiAnnotator(input_bedfile="../test_out/roc/roc_roc_psi_prediction.bed",
            annotation_bed6file="../test_data/annotation/hg38.genecode.v30.tRNA.snoRNA.miRNA.rmsk.exonFeatures.bed6",
            output_dir="../test_out/roc_annotation",
            output_name = "roc_annotation")

#annotate prediction result generated by SVM model
RpsiAnnotator(input_bedfile="../test_out/svm_predict/svm_predict_svm_psi_prediction.bed",
            annotation_bed6file="../test_data/annotation/hg38.genecode.v30.tRNA.snoRNA.miRNA.rmsk.exonFeatures.bed6",
            output_dir="../test_out/svm_annotation",
            output_name = "svm_annotation")

#annotate prediction result generated by ANN model
RpsiAnnotator(input_bedfile="../test_out/ann_predict/ann_predict_ann_psi_prediction.bed", 
            annotation_bed6file="../test_data/annotation/hg38.genecode.v30.tRNA.snoRNA.miRNA.rmsk.exonFeatures.bed6",
            output_dir="../test_out/ann_annotation",
            output_name = "ann_annotation")
```

#### Ψ-site target prediction
To run RpsiACAscan, we need to firstly prepare modification information data set.
```R
mod_info(RpsiScan_res_file="../test_out/user_defined/user_defined_ud_psi_prediction.bed",
         RpsiAnnotator_res_file="../test_out/ud_annotation/ud_annotation_anno.bed",
         psiU.SingleSites.bed="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
         pred_method="user-defined",
         output_dir="../test_out/mod_info",
         output_name="mod_info")

mod_info(RpsiScan_res_file="../test_out/roc/roc_roc_psi_prediction.bed",
         RpsiAnnotator_res_file="../test_out/roc_annotation/roc_annotation_anno.bed",
         psiU.SingleSites.bed="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
         pred_method="roc",
         output_dir="../test_out/mod_info",
         output_name="mod_info")

mod_info(RpsiScan_res_file="../test_out/svm_predict/svm_predict_svm_psi_prediction.bed",
         RpsiAnnotator_res_file="../test_out/svm_annotation/svm_annotation_anno.bed",
         psiU.SingleSites.bed="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
         pred_method="svm",
         output_dir="../test_out/mod_info",
         output_name="mod_info")

mod_info(RpsiScan_res_file="../test_out/ann_predict/ann_predict_ann_psi_prediction.bed",
         RpsiAnnotator_res_file="../test_out/ann_annotation/ann_annotation_anno.bed",
         psiU.SingleSites.bed="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
         pred_method="ann",
         output_dir="../test_out/mod_info",
         output_name="mod_info")
```

Secondly, we run psiACAscan to generate H/ACA snoRNA-RNA interaction information.

```R
#run psiACAscan to generate H/ACA snoRNA-RNA interaction information for ANN model prediction result
RpsiACAscan(acaboxseq="../test_data/annotation/human_hg38_snoRNABase_snoDB_rmRepeat.collapse.fa",
            modifiedfile="../test_out/mod_info/ann_predict_ann_psi_prediction_modinfo.txt",
            output_dir="../test_out/RpsiACAscan_out",
            output_name = "RpsiACAscan_out")

#organize psiACAscan result into tidy data table
tidy_psiACAscan(infile="../test_out/RpsiACAscan_out/RpsiACAscan_out.txt",
            appendfile="../test_out/mod_info/ann_annotation_anno.bed6.append",
            orphaninfoFile="../test_data/annotation/human_hg38_snoRNABase_snoDB_rmRepeat_addorphaninfo.csv",
            psiU.SingleSites.bed="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
            output_dir="../test_out/tidy_psiACAscan_out",
            output_name = "tidy_psiACAscan_out")
```

