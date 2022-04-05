##
library(rChoiceDialogs)
library(stringr)
select <- jchoose.files(default = "/media/respaldo4t/")
if(length(select)>0){
   bam.file <-RunSTARforARRIBA(sbjFile = select, nThreads = 10L)
   RunARRIBA(bam.file)
   sorted.bam.file <- RunSortIndexBam(bam.file, remove = F)
   FusionPlot(sorted.bam.file, savePlot = T)
 }else{
   cat("Not selected")
 }
