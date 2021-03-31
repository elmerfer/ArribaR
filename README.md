# ArribaR
ArribaR is an R based platform for Gene Expression and Gene Fusion based on the STAR aligner and the Arriba fusion detector.
It isolates all the requirements to run both softwares in an easy and transparent fashion, simplifying the instalation process and excusion.
It automatically download and install the basict softwares (STAr + Arriba) and all the required extra files like the GRCh37 assembly and GENCODE19 annotation as suggested by Arriba.

## Getting Started


## Installation
```
install.packages("devtools")
library(devtools)
install_github("elmerfer/ArribaR")
```

## Instaling software dependencies
```
InstallArribaR("where/I/want/to/instal/path")
```
This make take a long time and it will do the following process
* minimal installation of the STAR version 2.7.6a (the one required for Arriba version 2.1.0)
* the GRCc37(Hg19) Human reference file 
* the GENCODE GTF annotation file
* build the human reference index file
* minimal installation of the Arriba software 2.1.0
as a result it will leave an organized "Software" folder at where/I/want/to/instal/path

## Usage
```
runSTAR("my_1.fastq.gz") 
```

# Authors
