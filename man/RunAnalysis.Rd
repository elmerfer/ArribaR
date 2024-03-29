% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runners.R
\name{RunAnalysis}
\alias{RunAnalysis}
\title{RunAnalysis}
\usage{
RunAnalysis(sbjFile, feature, immuneTME = TRUE)
}
\arguments{
\item{sbjFile}{(character) full path to the first fasta file from a subject (subject_1.fastq or subject_1.fastq.gz ), the scond read file should be subject_2.fastq}

\item{feature}{(character)  "exon" or "gene". it will be overwritten with "gene" if immuneTME == TRUE}

\item{immuneTME}{(logical (default TRUE) should the \code{\link[MIXTURE]{MIXTURE}} algorithm be run?}
}
\value{
a list with the following slots:

Fusions : a Gene Fusion data (also in an excel file see \code{\link{RunArriba}})

Counts : a list as returned by \code{\link{GetCounts}}

iTME : immune tumor microenvironment objets as returned by \code{\link{GetImmuneContent}}
}
\description{
An all in one line to run STAR+Arriba+MIXTURE
}
\examples{
\dontrun{
results <- RunAnalysis(sbjFile = subject_1.fastq.gz)
##if a do not want to run TME and get the exon counts
results <- RunAnalysis(sbjFile = subject_1.fastq.gz, feature = "exon" , immuneTME = FALSE)
}
}
