\name{plot_psi_pos}
\alias{plot_psi_pos}
\title{Plot the pseudouridylation position}
\usage{
plot_psi_pos <- function(svgfile,
                         seqfile,
                         psifile,
                         output_dir,
                         output_name,
                         circle_radius="1",
                         circle_fill="red",
                         circle_stroke="black",
                         circle_stroke_width="0.1"
                         )
}
\arguments{
  \item{svgfile}{The svg file of RNA 2D structure (could be downloaded from RNAcentral).}
  \item{seqfile}{The sequence file of RNA (could be downloaded from RNAcentral).}
  \item{psifile}{The psi position file (pseudouridylation position of RNA, file contains extended sequence).}
  \item{circle_radius}{The radius of the circles.}
  \item{circle_fill}{The fill color of the circles (wrap the psi sites).}
  \item{circle_stroke}{The stroke of the circles (wrap the psi sites).}
  \item{circle_stroke_width}{The stroke width of the circles (wrap the psi sites).}
  \item{output_dir}{Directory of the output file.}
  \item{output_name}{Name of the output file.}
}
\description{This function will map the pseudouridylation sites in RNA 2D structures.}
\examples{

plot_psi_pos(svgfile="../test_data/plot_psi_pos_test/5_8S_rRNA_2D.svg",
             seqfile="../test_data/plot_psi_pos_test/5_8S_rRNA_sequence.txt",
             psifile="../test_data/plot_psi_pos_test/5_8S_rRNA_extended_sequencce.txt",
             output_dir="../test_out/5_8S_2D_rRNA",
             output_name="5_8S_2D_rRNA",
             circle_radius="1",
             circle_fill="red",
             circle_stroke="black",
             circle_stroke_width="0.1")

plot_psi_pos(svgfile="../test_data/plot_psi_pos_test/18S_rRNA_2D.svg",
             seqfile="../test_data/plot_psi_pos_test/18S_rRNA_sequence.txt",
             psifile="../test_data/plot_psi_pos_test/18S_rRNA_extended_sequencce.txt",
             output_dir="../test_out/18S_2D_rRNA",
             output_name="18S_2D_rRNA",
             circle_radius="1.5",
             circle_fill="red",
             circle_stroke="black",
             circle_stroke_width="0.1")

plot_psi_pos(svgfile="../test_data/plot_psi_pos_test/28S_rRNA_2D.svg",
             seqfile="../test_data/plot_psi_pos_test/28S_rRNA_sequence.txt",
             psifile="../test_data/plot_psi_pos_test/28S_rRNA_extended_sequencce.txt",
             output_dir="../test_out/28S_2D_rRNA",
             output_name="28S_2D_rRNA",
             circle_radius="1.1",
             circle_fill="red",
             circle_stroke="black",
             circle_stroke_width="0.1")

}
\value{Generate figure of RNA 2D structures with pseudouridylation sites highlight.}
