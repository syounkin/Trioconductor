\name{FamilyExperiment-class}
\docType{class}

% Class:
\alias{class:FamilyExperiment}
\alias{FamilyExperiment-class}
\alias{FamilyExperiment}

% Accessors:
\alias{geno,FamilyExperiment-method}

% FE methods:
\alias{MAF,FamilyExperiment-method}

% "show" method:
\alias{show,FamilyExperiment-method}

\title{FamilyExperiment Class}

\description{The FamilyExperiment (FE) class is an extension of the
  SummarizedExperiment (SE) class defined in GenomicRanges. An FE object
  contains an additional slot for a pedigree object.  If the FE object
  is designed to contain genotype data, then a SNPMatrix object named
  'geno' must appear in the assayData list.  Similarly with cnv data the
  assayData list must contain a matrix named 'cnv'.}

\author{Samuel G. Younkin <syounkin@jhsph.edu>}
