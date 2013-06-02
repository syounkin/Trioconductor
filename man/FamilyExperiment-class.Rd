\name{FamilyExperiment-class}
\docType{class}

% Class:
\alias{class:FamilyExperiment}
\alias{FamilyExperiment-class}
\alias{FamilyExperiment}

% Accessors:
\alias{pedigree,FamilyExperiment-method}
\alias{pedigree}
\alias{MAF,FamilyExperiment-method}
\alias{MAF}
\alias{completeTrios,FamilyExperiment-method}
\alias{completeTrios}
\alias{parents,FamilyExperiment-method}
\alias{parents}

\title{FamilyExperiment Class}

\description{The FamilyExperiment (FE) class is an extension of the
  SummarizedExperiment (SE) class defined in GenomicRanges. An FE object
  contains an additional slot for a pedigree object.  If the FE object
  is designed to contain genotype data, then a SNPMatrix object named
  'geno' must appear in the assayData list.  Similarly with cnv data the
  assayData list must contain a matrix named 'cnv'.}

  \section{Accessors}{
  In the following code snippets, \code{fe} is a FamilyExperiment
  object.
  \describe{
    \item{}{
      \code{pedigree(fe)}: Returns the pedigree object
    contained in \code{fe}.  The returned object will have class
    PedClass.
  }
  \item{}{
    \code{MAF(fe)}: Returns a numeric vector of minor allele
    frequencies.  The frequencies are computed after subsetting the data
  to be parents only.  This clearly is problematice if your pedigree
  contains multiple generations.  This method is essentially a wrapper
  for the The minor allele frequency method of SnpStats.
}
\item{}{
\code{completeTrios(fe)}: Returns an \code{n} by 3 (character or factor)?
  matrix with rows that correspond to the parent-offspring trios
  contained in \code{fe} such that each of the trio members appears in
  the assay data.  \code{n} is the number of ``complete trios.''
  }
\item{}{
\code{parents(fe)}: Returns the parents IDs.  Note that this is
  intended to be used with parent-offspring trios.  It will return any
  subject that has both parents in the pedigree data.
  }
}
}
\author{Samuel G. Younkin <syounkin@jhsph.edu>}
