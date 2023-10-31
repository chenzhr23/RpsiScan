\name{ground_truth}
\alias{ground_truth}
\title{This is a function to prepare ground truth data set for pseudouridylation identification}
\usage{
ground_truth <- function(RpsiScan_res_file,
                         known_rRNA_bed,
                         rrna_chr21_bed,
                         output_dir,
                         output_name)
}
\arguments{
  \item{RpsiScan_res_file}{The file path of psiFidner result (RpsiScan output file).}
  \item{known_rRNA_bed}{The bed file contains knonw rRNA pseudouridylation sites, default we recommend hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed.}
  \item{rrna_chr21_bed}{The bed file contains rRNA region in chromosome 21 (18S/28S/5.8S), default we recommend rrna_chr21.bed.}
  \item{output_dir}{Directory of the output file.}
  \item{output_name}{Name of the output file.}
}
\description{
Prepare ground truth data set (known rRNA positive/negative pseudouridylation sites) for pseudouridylation identification.
}
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
}
\value{
tab-delimited text files include ground truth rRNA (based on known rRNA pseudouridylation sites) reverse transcription stop read information for each CMC-treatment and CMC-control group, e.g. treatStopRatio, ctrlStopRatio
}
