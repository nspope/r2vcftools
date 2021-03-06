% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion.R
\docType{methods}
\name{GenotypeMatrix}
\alias{GenotypeMatrix}
\alias{GenotypeMatrix,vcfLink-method}
\title{Convert VCF to genotype matrix and return}
\usage{
GenotypeMatrix(object, ...)

\S4method{GenotypeMatrix}{vcfLink}(object, filename = NULL, verbose = TRUE)
}
\arguments{
\item{object}{an object of class "vcfLink"}

\item{verbose}{print stdout from vcftools?}

\item{output.file}{the filename to save the genotype matrix to, if you want to write it outside of R}
}
\description{
Given a "vcfLink" object, converts the VCF file to genotype matrix format and returns this matrix. See the "--012" flag in the vcftools documentation.
}
