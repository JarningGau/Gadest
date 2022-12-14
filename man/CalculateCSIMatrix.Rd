% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summerize_DEGs.R
\name{CalculateCSIMatrix}
\alias{CalculateCSIMatrix}
\title{Calculating connection specificity index (CSI) matrix.}
\usage{
CalculateCSIMatrix(pccMat, cores = 1)
}
\arguments{
\item{pccMat}{pearson correlation coefficient (PCC) matrix of genes.}

\item{cores}{the number of threads used. Default: 1}
}
\value{
a CSI matrix
}
\description{
Calculating connection specificity index (CSI) matrix.
}
