% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LTMRmain.R
\name{LTMmeta}
\alias{LTMmeta}
\title{Integrate the LTM analysis results of multiple datasets.}
\usage{
LTMmeta(metaL, fileName = NULL, outDir = NULL)
}
\arguments{
\item{metaL}{a list of data frame. Each data frame is from LTMdish output.}

\item{fileName}{a character string. The name of output file. Set it \code{NULL}(default) if prefer to return a data frame.}

\item{outDir}{a character string. The name of directory used to store the output file.}
}
\value{
a data fame
}
\description{
This function get the meta values by integrating LTM analysis results of multiple datasets, e.g., multiple skin datasets.
}
\examples{
## please refer to the webpage of LTMR package
}
