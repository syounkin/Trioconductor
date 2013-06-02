\name{PedClass-class}
\docType{class}

% Class:
\alias{class:PedClass}
\alias{PedClass-class}
\alias{PedClass}

% Accessors:
\alias{offspring,PedClass-method}
\alias{offspring}
\alias{parents,PedClass-method}
\alias{parents}
\alias{trios,PedClass-method}
\alias{trios}
\alias{allSubjects,PedClass-method}
\alias{allSubjects}

\title{PedClass Class}

\description{Container for pedigree data}

  \section{Accessors}{
  In the following code snippets, \code{ped} is a PedClass
  object.
  \describe{
    \item{}{
      \code{offspring(ped)}: Returns a (character or factor)? vector of
  offspring ID.
  }
    \item{}{
      \code{parents(ped)}: Returns a (character or factor)? vector of
  parents ID.
  }
    \item{}{
      \code{trios(ped)}: Returns a (character or factor)? vector of
  trios ID.
}
\item{}{
\code{allSubjects(fe)}: Returns a (character or factor)? vector
  of the subject IDs contained in the pedigree object of
  \code{fe}
  }
}
}
\author{Samuel G. Younkin <syounkin@jhsph.edu>}
