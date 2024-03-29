% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runners.R
\name{RunSTAR}
\alias{RunSTAR}
\title{RunSTAR
This function will run STAR software}
\usage{
RunSTAR(sbjFile, nThreads, assemblyVersion = c("37", "38"))
}
\arguments{
\item{sbjFile}{string with the full path and name of the xxxx_1_fastq or gz sequence file. This
function only support paired data}

\item{nThreads}{number of CPUs threads}
}
\description{
RunSTAR
This function will run STAR software
}
\examples{
\dontrun{
test.subject <- GetArribaRTest()
bam.subject <- RunSTAR(test.subject)
}
}
\author{
Elmer A. Fernández
}
