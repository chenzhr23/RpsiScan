\name{prefilt}
\alias{prefilt}
\title{Filt out undesired psiScan information}
\usage{
prefilt <- function(RpsiScan_res,
                    target_chrom="^chr[0-9|a-z|A-Z]*$",
                    target_base="T",
                    output_dir,
                    output_name)
}
\arguments{
  \item{RpsiScan_res}{The R result of RpsiScan function.}
  \item{target_chrom}{The target chromsome to be retained, default is "^chr[0-9|a-z|A-Z]*$", meaning retain those psiScan information with normal chromsome.}
  \item{target_base}{The target base to be retained, default is "T", meaning retain those psiScan information with T base (where pseudouridylation take place).}
  \item{output_dir}{Directory of the output file.}
  \item{output_name}{Name of the output file.}
}
\description{This function will retain psiScan information with desired target chrom/target base.}
\examples{
RpsiScan_example()#when finished, view your global environment, and you will detect Data named 'RpsiScan_example'

prefilt_out <- prefilt(RpsiScan_res=RpsiScan_example,
                        target_chrom="^chr[0-9|a-z|A-Z]*$",
                        target_base="T",
                        output_dir="../test_out/RpsiScan_out",
                        output_name="test_prefilt")

}
\value{Retained desired psiScan information (by pre-filt argument)}
