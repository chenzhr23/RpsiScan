# RpsiScan_example
#
# Get RpsiScan example result

#' Get psiScan example result
#' @export
#'
RpsiScan_example <- function()
{
  #1. /public/home/chenzr/miniconda3/lib/R/library/RpsiScan/extdata/psiScan --fa ./test_data/genome/hg38.fa --fai ./test_data/genome/hg38.fa.fai --treat ./test_data/A2.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gzAligned.sortedByCoord.out.bam --input ./test_data/A1.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gzAligned.sortedByCoord.out.bam -p 1.5 -t 5 -r 0.05 -M 1 -f 1 -m 0 -s -n -w 20 --gene ./test_data/annotation/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12 -o ./test_out/test.txt 2>./test_out/test.txt.log
  #2. psiScan_res_filt<-psiScan_res %>% filter(base=="T",treatStopNum>10)
  #3. fwrite(psiScan_res_filt,"CIV-seq_test.txt")
  assign("RpsiScan_example",fread(system.file("example", "CIV-seq_test.txt", package = "RpsiScan")),envir = .GlobalEnv)
  print("This is a RpsiScan example result generated from CIV-seq samples.")

}
