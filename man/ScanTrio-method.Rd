\name{ScanTrio-method}
\docType{methods}

\alias{ScanTrio}
\alias{ScanTrio,FamilyExperiment,GRanges,GRanges-method}

\title{ScanTrio method}

\description{This method implements a trio-based version of the Scan
  Statistic approach of Ionita-Laza et al.}

\usage{
\S4method{ScanTrio}{FamilyExperiment,GRanges,GRanges}(object, window, block)
}

\arguments{
  \item{object}{
    A \link{FamilyExperiment} object. 
  }
  \item{window}{
    A GRanges object. 
  }
  \item{block}{
    A GRanges object. 
  }
}
\details{Foo.}
\value{A GRanges object.}
\author{Samuel G. Younkin <syounkin@jhsph.edu>}
