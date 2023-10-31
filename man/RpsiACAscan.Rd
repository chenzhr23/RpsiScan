\name{RpsiACAscan}
\alias{RpsiACAscan}
\title{R API for psiACAscan}
\usage{
RpsiACAscan <- function(acaboxseq,
                        modifiedfile,
                       output_dir,
                       output_name)
}
\arguments{
  \item{acaboxseq}{Sequence file in fasta format, default we recommend human_hg38_snoRNABase_snoDB_rmRepeat.collapse.fa.}
  \item{modifiedfile}{Modified file with three column header: target_site_name (e.g. chr1_10299004_10299005_+), snoRNA_information (e.g. U6-40-HBI-100/ACA12), target_site_sequence (e.g. CUAAAAUUGGAACGAYACAGAGAAGAUUAGC).}
  \item{output_dir}{Directory of the output file.}
  \item{output_name}{Name of the output file.}
  }
\description{RpsiACAscan will run psiACAscan to generate H/ACA snoRNA-RNA interaction information.}
\examples{
RpsiACAscan(acaboxseq="../test_data/annotation/human_hg38_snoRNABase_snoDB_rmRepeat.collapse.fa",
            modifiedfile="../test_out/mod_info/ann_predict_ann_psi_prediction_modinfo.txt",
            output_dir="../test_out/RpsiACAscan_out",
            output_name = "RpsiACAscan_out.txt")
}
\value{Interaction information between H/ACA snoRNAs and target RNAs}

