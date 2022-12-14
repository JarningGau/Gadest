% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_expression_shift.R
\name{sample.nmi}
\alias{sample.nmi}
\title{Calculating median NMI for a given cells.}
\usage{
sample.nmi(
  X,
  cells.1,
  cells.2 = NULL,
  N.cells = 100,
  seed = 1,
  sample.on.cell.pairs = TRUE,
  N.cell.pairs = 100
)
}
\arguments{
\item{X}{matrix, discrete gene expression matrix. Rows, genes; columns, cells.}

\item{cells.1}{vector, selected cells for group1. Match to colnames(X)}

\item{cells.2}{vector, selected cells for group2. Match to colnames(X)}

\item{N.cells}{int, sample size of selected cells.}

\item{seed}{int, seed for sampling.}

\item{sample.on.cell.pairs}{bool, whether sample on cell pairs for NMI calculation.}

\item{N.cell.pairs}{int, sample size of cell pairs for NMI calculation.}
}
\value{
numeric, median NMI
}
\description{
Calculating median NMI for a given cells.
}
