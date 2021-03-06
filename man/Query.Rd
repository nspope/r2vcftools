% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\docType{methods}
\name{Query}
\alias{Query}
\alias{Query,vcfLink-method}
\title{Query attributes of a VCF file}
\usage{
Query(object, ...)

\S4method{Query}{vcfLink}(object, type = c("counts2", "freq2", "site-depth",
  "site-mean-depth", "geno-depth", "depth", "site-quality", "het", "hardy",
  "site-pi"), verbose = TRUE)
}
\arguments{
\item{object}{an object of class "vcfLink"}

\item{type}{the type of information to return; see the vcftools manual for details}

\item{verbose}{boolean; if true, return stdout from vcftools}
}
\description{
Use vcftools to extract information from a VCF file.
}
