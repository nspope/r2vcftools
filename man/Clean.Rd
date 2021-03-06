% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{methods}
\name{Clean}
\alias{Clean}
\alias{Clean,vcfLink-method}
\title{Clean unused temporary VCF files}
\usage{
Clean(object, ...)

\S4method{Clean}{vcfLink}(object)
}
\arguments{
\item{object}{an object of class "vcfLink"}
}
\description{
Sweeps the temporary directory and clears temporary VCF files that are not currently linked to an object.
}
