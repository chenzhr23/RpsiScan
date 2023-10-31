# plot_psi_pos
#
# Plot the pseudouridylation position

#' Plot the pseudouridylation position
#' @param svgfile The svg file of RNA 2D structure (could be downloaded from RNAcentral)
#' @param seqfile The sequence file of RNA (could be downloaded from RNAcentral)
#' @param psifile The psi position file (pseudouridylation position of RNA, file contains extended sequence)
#' @param circle_radius The radius of the circles
#' @param circle_fill The fill color of the circles (wrap the psi sites)
#' @param circle_stroke The stroke of the circles (wrap the psi sites)
#' @param circle_stroke_width The stroke width of the circles (wrap the psi sites)
#' @param output_dir The path to the output directory
#' @param output_name The output file name
#' @export
#'
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
{
  # svgfile="28S_rRNA_2D.svg"
  # seqfile="28S_rRNA_sequence.txt"
  # psifile="28S_rRNA_extended_sequencce.txt"
  # output_dir="./28S_2D_rRNA"
  # output_name="28S_2D_rRNA"

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  output_dir_output_name<-paste(output_dir,"/",output_name,sep="")

  ####Get psi position####

  seqfile <- file(seqfile, open = "r")
  seqfile_in <- readLines(seqfile)
  close(seqfile)

  psifile <- file(psifile, open = "r")
  psifile_in <- readLines(psifile)
  close(psifile)

  tmp_pos<-vector()
  lapply(seq_along(psifile_in),function(i){
    tmp_pos[i]<<- gregexpr(psifile_in[i], seqfile_in)[[1]]+15
  })

  ####Modify svg file####
  svgfile_con <- file(svgfile, open = "r")
  RNA_2D_structure <- readLines(svgfile_con)
  close(svgfile_con)
  RNA_2D_structure1<-sub("<g><style type=\"text/css\" >","<g><style type=\"text/css\" >  .stA{fill:#258B3A;} .stG{fill:#C26E19;} .stC{fill:#1A4693;} .stU{fill:#D9291D;}",RNA_2D_structure)
  RNA_2D_structure2<-sub("class=\"green\" >5\'</text></g>","class=\"numbering-label sequential\" >5\'</text></g>",RNA_2D_structure1)
  RNA_2D_structure3<-sub("class=\"green\" >5\'</text></g>","class=\"numbering-label sequential\" >5\'</text></g>",RNA_2D_structure2)
  RNA_2D_structure4<-sub("class=\"green\" >3\'</text></g>","class=\"numbering-label sequential\" >3\'</text></g>",RNA_2D_structure3)
  RNA_2D_structure5<-sub("class=\".*\" >C</text></g>","class=\"stC\" >C</text></g>",RNA_2D_structure4)
  RNA_2D_structure6<-sub("class=\".*\" >A</text></g>","class=\"stA\" >A</text></g>",RNA_2D_structure5)
  RNA_2D_structure7<-sub("class=\".*\" >G</text></g>","class=\"stG\" >G</text></g>",RNA_2D_structure6)
  RNA_2D_structure8<-sub("class=\".*\" >U</text></g>","class=\"stU\" >U</text></g>",RNA_2D_structure7)

  svgfile_coordinate<-read.delim(svgfile,header = F)
  svgfile_coordinate_V1<-data.frame(row_index=svgfile_coordinate$V1)
  svgfile_coordinate_df<-svgfile_coordinate_V1 %>% filter(grepl("<g><title>",row_index))
  svgfile_coordinate_df_sep<-svgfile_coordinate_df %>% separate_wider_delim(row_index,delim ="=" ,names=c("tag1","tag2","tag3","tag4"))
  svgfile_coordinate_df_sep$x_pos<-sub(" y","",svgfile_coordinate_df_sep$tag2)
  svgfile_coordinate_df_sep$y_pos<-sub(" class","",svgfile_coordinate_df_sep$tag3)
  svgfile_coordinate_df_sep$base<-sub("</text></g>","",svgfile_coordinate_df_sep$tag4)
  svgfile_coordinate_df_sep$base<-sub(".*>","",svgfile_coordinate_df_sep$base)
  svgfile_coordinate_df_sep$pos<-sub(" (.*)</title><text x","",svgfile_coordinate_df_sep$tag1)
  svgfile_coordinate_df_sep$pos<-sub("<g><title>","",svgfile_coordinate_df_sep$pos)
  svgfile_coordinate_df_sep<-svgfile_coordinate_df_sep %>% filter(!grepl("numbering-label sequential.*",svgfile_coordinate_df_sep$tag4))
  svgfile_coordinate_df_sep_tmp_pos<-svgfile_coordinate_df_sep %>% filter(grepl(paste(paste("^",tmp_pos, "$", sep= ""),collapse = "|"),svgfile_coordinate_df_sep$pos))
  # <circle cx="81.3899" cy="98.0298" r="1" fill="red" stroke="red" stroke-width="0.1"/>
  add_circle_text<-paste("<circle cx=",
                         "\"",
                         svgfile_coordinate_df_sep_tmp_pos$x_pos,
                         "\"",
                         " cy=",
                         "\"",
                         svgfile_coordinate_df_sep_tmp_pos$y_pos,
                         "\"",
                         " r=",
                         "\"",
                         circle_radius,
                         "\"",
                         " fill=",
                         "\"",
                         circle_fill,
                         "\"",
                         " stroke=",
                         "\"",
                         circle_stroke,
                         "\"",
                         " stroke-width=",
                         "\"",
                         circle_stroke_width,
                         "\"",
                         " />",
                         sep="")

  RNA_2D_structure9<-c(RNA_2D_structure8[-length(RNA_2D_structure8)],paste(paste(add_circle_text,collapse="\n"),"</svg>",sep=""))

  svgout <- file(paste(output_dir_output_name,".svg",sep=""), open = "w")
  writeLines(RNA_2D_structure9, svgout)
  close(svgout)
}
