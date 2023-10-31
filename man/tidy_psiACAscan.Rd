\name{tidy_psiACAscan}
\alias{tidy_psiACAscan}
\title{Organize psiACAscan result into tidy data table}
\usage{
tidy_psiACAscan <- function(infile,
                            appendfile,
                            orphaninfoFile,
                            psiU.SingleSites.bed,
                            output_dir,
                            output_name)
}
\arguments{
  \item{inFile}{The input file generated from mod_info() function.}
  \item{appendfile}{The path to bed6.append file generated from mod_info() function.}
  \item{orphaninfoFile}{The path to known snoRNA information file.}
  \item{psiU.SingleSites.bed}{psiU.SingleSites.bed pseudouridylation sites in bed format, default we recommend hg38.psiU.SingleSites.bed.}
  \item{output_dir}{Directory of the output file.}
  \item{output_name}{Name of the output file.}
  }
\description{tidy_psiACAscan will run psiACAscan to generate H/ACA snoRNA-RNA interaction information.}
\examples{
mod_info(RpsiScan_res_file="../test_out/ann_predict/ann_predict_ann_psi_prediction.bed",
         RpsiAnnotator_res_file="../test_out/ann_annotation/ann_annotation_anno.bed",
         psiU.SingleSites.bed="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
         pred_method="ann",
         output_dir="../test_out/mod_info",
         output_name="mod_info")

RpsiACAscan(acaboxseq="../test_data/annotation/human_hg38_snoRNABase_snoDB_rmRepeat.collapse.fa",
            modifiedfile="../test_out/mod_info/ann_predict_ann_psi_prediction_modinfo.txt",
            output_dir="../test_out/RpsiACAscan_out",
            output_name = "RpsiACAscan_out.txt")

tidy_psiACAscan(infile="../test_out/RpsiACAscan_out/RpsiACAscan_out.txt",
            appendfile="../test_out/mod_info/ann_annotation_anno.bed6.append",
            orphaninfoFile="../test_data/annotation/human_hg38_snoRNABase_snoDB_rmRepeat_addorphaninfo.csv",
            psiU.SingleSites.bed="../test_data/ground_truth/hg38.psiU.SingleSites.bed",
            output_dir="../test_out/tidy_psiACAscan_out",
            output_name = "tidy_psiACAscan_out")
}
\value{Well-organized data tale of psiACAscan information contain H/ACA snoRNA-target RNA interaction}

