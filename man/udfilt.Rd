\name{udfilt}
\alias{udfilt}
\title{Filt out undesired psiScan information}
\usage{
udfilt <- function(RpsiScan_res_file,
                   tretRtsRatio_thres,
                   rtsRatioFold_thres,
                   tretBefRatio_thres,
                   befRatioFold_thres,
                   tretAftRatio_thres,
                   aftRatioFold_thres,
                   tretMutRatio_thres,
                   mutRatioFold_thres,
                   tretDelRatio_thres,
                   delRatioFold_thres,
                   output_dir,
                   output_name)
}
\arguments{
  \item{RpsiScan_res_file}{The file path of psiFidner result (RpsiScan output file).}
  \item{tretRtsRatio_thres}{The filtering threshold for tretRtsRatio.}
  \item{rtsRatioFold_thres}{The filtering threshold for rtsRatioFold.}
  \item{tretBefRatio_thres}{The filtering threshold for tretBefRatio.}
  \item{befRatioFold_thres}{The filtering threshold for befRatioFold.}
  \item{tretAftRatio_thres}{The filtering threshold for tretAftRatio.}
  \item{aftRatioFold_thres}{The filtering threshold for aftRatioFold.}
  \item{tretMutRatio_thres}{The filtering threshold for tretMutRatio.}
  \item{mutRatioFold_thres}{The filtering threshold for mutRatioFold.}
  \item{tretDelRatio_thres}{The filtering threshold for tretDelRatio.}
  \item{delRatioFold_thres}{The filtering threshold for delRatioFold.}
  \item{output_dir}{Directory of the output file.}
  \item{output_file}{Name of the output file.}
}
\description{This function will retain psiScan information with desired identification metrics greater than thresholds (treatPreRpmFold, treatAfterRpmFold, treatStopRatio, preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC).}
\examples{
RpsiScan_example()#when finished, view your global environment, and you will detect Data named 'RpsiScan_example'

prefilt_out <- prefilt(RpsiScan_res=RpsiScan_example,
                        target_chrom="^chr[0-9|a-z|A-Z]*$",
                        target_base="T",
                        output_dir="../test_out/RpsiScan_out",
                        output_name="test_prefilt")

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
}
\value{Retained desired psiScan information (by user-defined argument)}
