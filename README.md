# ArribaR
ArribaR is an R based platform for Gene Expression and Gene Fusion identification based on the [STAR](https://github.com/alexdobin/STAR) aligner and the [Arriba](https://arriba.readthedocs.io/en/latest/) fusion detector.
It is part of an analytical suite for the study of genomic rearrangements, its potential as target therapy and as neoantigenic role. In addition it may provide the immune microenvironment by [MIXTURE](https://github.com/elmerfer/MIXTURE) by the LM22 molecular signature
The ArribaR software package is intended to provide a friendly and effortless R tool for non bioinformatics. It isolates all the requirements to run both softwares in an easy and transparent fashion, simplifying the installation process and execution.
It automatically download and install the basic softwares (STAR + Arriba) extra files, like the genome references and annotations should be downloaded using the [GenomeDB](https://github.com/elmerfer/GenomeDB) and [Aligners](https://github.com/elmerfer/Aligners) R packages of this family Next Generation Sequencing platform.

## Getting Started
Please visit
[GenomeDB](https://github.com/elmerfer/GenomeDB) and [Aligners](https://github.com/elmerfer/Aligners).

*[Install the Human Genome ](https://github.com/elmerfer/GenomeDB/wiki)

*[Install and build the STAR reference index](https://github.com/elmerfer/Aligners/wiki)

## Install ArribaR package
```
install.packages("devtools")
library(devtools)
install_github("elmerfer/ArribaR")
```

## Install Arriba software dependencies
```
InstallArribaR()
```
<<<<<<< HEAD
The instalation will be performed under the [GenomeDB](https://github.com/elmerfer/GenomeDB) organization
=======
This may take a long time and it will do the following process
* minimal installation of the STAR version 2.7.6a (the one required for Arriba version 2.1.0)
* the GRCc37(Hg19) Human reference file 
* the GENCODE GTF annotation file
* build the human reference index file
* minimal installation of the Arriba software 2.1.0

as a result it will leave an organized "Software" folder at where/I/want/to/install/path/Software


## Usage
This implementation is specifically provided to process human RNAseq paired data, thus requiring to hold both xxxx_1.fastq(.gz) and xxxx_2.fastq(.gz).
The row sequencing data should be "xxxxx_1.fastq" or "xxxxx_1.fastq.gz". The second sequen file will be automatically defined as "xxxxx_2.fastq" or "xxxxx_2.fastq.gz"
```
library(ArribaR)
test.subject <- GetArribaRTest()
bam.subject <- RunSTARforARRIBA(test.subject)
Fusions <- RunARRIBA(bam.subject)
View(Fusions)
```
## Interactive usage
copy this R script in your RStudio and source it
```
library(stringr)
select <- file.choose()
if(length(select)>0){
   bam.file <-RunSTARforARRIBA(sbjFile = select, nThreads = 10L)
   RunARRIBA(bam.file)
   sorted.bam.file <- RunSortIndexBam(bam.file, remove = F)
   FusionPlot(sorted.bam.file, savePlot = T)
 }else{
   cat("No File selected")
 }
 ```

# Authors
* **Elmer A. Fernández** CIDIE-UCC-CONICET
