##
library(rChoiceDialogs)
library(ArribaR)
library(stringr)
select <- rstudioapi::selectDirectory(
  caption = "Select Directory",
  label = "Select",
  path = "/home/biomolecular/DATA/NGS/RNAseq/Pacientes"
)

if(length(select)>0){
  
  
   dir.create(file.path("/mnt/data/RNAseq/Patients",basename(select)))
   lfo <- list.files(select, full.names = T, pattern = ".fastq.gz")
   ok<-file.copy(from = lfo, to = file.path("/mnt/data/RNAseq/Patients",basename(select)))
   print(lfo)
   subjPath <- file.path("/mnt/data/RNAseq/Patients",basename(select))
  # subjPath<-select
  lf <- list.files(subjPath, full.names = T, pattern = "1.fastq")
  
   bam.file <-RunSTARforARRIBA(sbjFile = lf, nThreads = 10L)
   ArribaR::RunARRIBA(bam.file, allProteinPredictions = T)
   sorted.bam.file <- RunSortIndexBam(bam.file, remove = F)
   FusionPlot(sorted.bam.file, savePlot = T)
  print(lf)
  
  
 }else{
   cat("Not selected")
 }
 # FusionPlot(sorted.bam.file, savePlot = T)

