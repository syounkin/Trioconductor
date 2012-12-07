\name{TrioSet}
\alias{TrioSet}

\title{
  Constructor for TrioSet objects.
}
\description{
  Construct an object of class \code{TrioSet}.
}
\usage{
TrioSet(pedigreeData = Pedigree(), sample.sheet, row.names = NULL, lrr,
  baf, featureData, cdfname, drop = TRUE, mindist=NULL, genome=c("hg19",
  "hg18"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{pedigreeData}{
    A \code{Pedigree}  object.
  }
  \item{sample.sheet}{
    A \code{data.frame} of sample covariates.
  }
  \item{row.names}{

      Row identifiers for \code{sample.sheet} that match the names of
   the
   trios in the \code{Pedigree} object.

 }
 \item{lrr}{
   Matrix of log R ratios of dimension p by 3n where n is the numper of
 offspring and p is the number of markers.  The order of the columns
 must be F, M, O for each trio.  To create such a matrix from the
 pedigree and named vector of log R ratios see trioMat().
 }

 \item{baf}{
   Matrix of B allele frequencies of dimension p by 3n where n is the numper of
 offspring and p is the number of markers.  The order of the columns
 must be F, M, O for each trio.  To create such a matrix from the
 pedigree and named vector of log R ratios see trioMat().
 }

 \item{featureData}{
   \code{AnnotatedDataFrame} for features.
 }

 \item{cdfname}{
   \code{character} string for annotation package.
 }

 \item{drop}{
       \code{Logical}. If FALSE, the row and column \code{dimnames} are
  null.
}

\item{mindist}{
  A \code{matrix} of the minimum distance.
}

\item{genome}{\code{character}: the UCSC genome build used for the
  probe annotation}

}
