# ArribaR
ArribaR is an R based platform for Gene Expression and Gene Fusion identification based on the STAR aligner and the Arriba fusion detector.
It is part of an analitycal suite for the study of genomic rearrangements, its potential as target therapy, their neoantigenic roles and the immune microenvironment by [MIXTURE](https://github.com/elmerfer/MIXTURE)
The ArribaR software package is intended to provide a friendly and effortless R tool for non bioinformatics. It isolates all the requirements to run both softwares in an easy and transparent fashion, simplifying the installation process and execusion.
It automatically download and install the basic softwares (STAr + Arriba) and all the required extra files like the GRCh37 assembly and GENCODE19 annotation as suggested by Arriba.

## Getting Started


## Installation
```
install.packages("devtools")
library(devtools)
install_github("elmerfer/ArribaR")
```

## Instaling software dependencies
```
InstallArribaR("where/I/want/to/install/path")
```
This make take a long time and it will do the following process
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
runAnalisis("my_1.fastq.gz") 
```

# Authors
* **Elmer A. Fernández** CIDIE-UCC-CONICET
