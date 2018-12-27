markDown_in=$1
pdf_Out=$2
pandoc_path="/nas/longleaf/apps/rstudio/1.0.136/bin/pandoc"


R -e Sys.setenv"(RSTUDIO_PANDOC='$pandoc_path')" -e  rmarkdown::render"('$markDown_in',output_file='$pdf_Out')"