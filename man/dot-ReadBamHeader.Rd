% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runners2.R
\name{.ReadBamHeader}
\alias{.ReadBamHeader}
\title{.ReadBamHeader}
\usage{
.ReadBamHeader(bamFile)
}
\arguments{
\item{bamFile}{bam file full path}
}
\value{
a list with the following slots:
CHR: the list of chromosomes and sequence names in the bam file
ProgramName : the alignment program (PG:PN see [BAM header](https://samtools.github.io/hts-specs/SAMv1.pdf))
ProgramVersion : the version of the aligner
Code : executed source line code
GenomeDBversion : the genome version used (see GenomeDB)
GenomeDBpath : the path to genome version
}
\description{
.ReadBamHeader
}
