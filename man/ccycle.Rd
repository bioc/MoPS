\name{ccycle}
\alias{ccycle}
\docType{data}
\title{
Dataset containing 500 time series gene expression measurements.
}
\description{
Example dataset as used in the MoPS Vignette Case Study.
Data matrix with 500 time series of gene expression measurements of 41 consecutive measurements seperated by 5 minutes (columns).
It is a subset of the dataset published in Eser et al. (Mol Sys Biol, 2014). ArrayExpress accession: E-MTAB-1908.     
}
\usage{data(ccycle)}
\format{
  The format is:
 num [1:500, 1:41] 
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:500] "YKR077W" "YJL218W" "YGR009C" "YIL104C" ...
  ..$ : chr [1:41] "5" "10" "15" "20" ...
}
\details{
unique gene identifiers are given as row names. 
}
\source{
Periodic mRNA synthesis and degradation co-operate during cell cycle gene expression (Eser et al. Mol Sys Biol, 2014) 
}
\keyword{datasets}
