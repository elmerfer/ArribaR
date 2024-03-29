% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runners2.R
\name{RunSTARforARRIBA}
\alias{RunSTARforARRIBA}
\title{RunSTARforARRIBA
This function will run STAR software
based on GenomeDB}
\usage{
RunSTARforARRIBA(sbjFile, nThreads, version = "GRCh38+GENECODE")
}
\arguments{
\item{sbjFile}{string with the full path and name of the xxxx_1_fastq or gz sequence file. This
function only support paired data}

\item{nThreads}{number of CPUs threads}
}
\description{
RunSTARforARRIBA
This function will run STAR software
based on GenomeDB
}
\details{
In order to look for gene fusions with the ARRIBA software, the STAR aligner should be run with special parameters.
The RunSTARforARRIBA function performe such task.
}
\examples{
\dontrun{
test.subject <- GetArribaRTest()
bam.subject <- RunSTAR2(test.subject)
##to display the genome and annotation version
attr(out.file,"GenomeDB")
##to display the elapsed time of the STAR run
attr(out.file,"ElapsedTime") 
}
}
\author{
Elmer A. Fernández
}
