\name{psiSVMpred}
\alias{psiSVMpred}
\title{Predict pseudouridylation sites based on SVM model}
\usage{
psiSVMpred <- function(filtfile,
                       svmmodelfile,
                       rRNAfile,
                       output_dir,
                       output_name)
}
\arguments{
  \item{filtfile}{File to be filted (File ready to be predicted).}
  \item{svmmodelfile}{SVM model file.}
  \item{rRNAfile}{The bed file contains knonw rRNA pseudouridylation sites, default we recommend hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed.}
  \item{output_dir}{Directory of the output file.}
  \item{output_name}{Name of the output file.}
}
\description{This function will predict pseudouridylation sites based on SVM model.}
\examples{
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

psiSVM(rocfile="../test_out/ground_truth/ground_truth_roc_plot.txt",
       rRNAfile="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
       rRNAfile2="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
       filtfile="../test_out/RpsiScan_out/test_prefilt.txt",
       output_dir="../test_out/svm_model",
       output_name="svm_model")

psiSVMpred(filtfile="../test_out/RpsiScan_out/test_prefilt.txt",
           svmmodelfile="../test_out/svm_model/svm_model_SVM_model.RData",
           rRNAfile="../test_data/ground_truth/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed",
           output_dir="../test_out/svm_predict",
           output_name="svm_predict")

}
\value{Predicted pseudouridylation result based on SVM model}
