% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dna.R
\name{dna_pathway}
\alias{dna_pathway}
\title{Differential network analysis on a single pathway}
\usage{
dna_pathway(
  x,
  pathway,
  groups = NULL,
  network_inference = run_pcor,
  n_perm = 100,
  lp = 2,
  seed = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{The gene expression data to be analyzed. Assumed to already be
processed by \code{\link{dnapath}}}

\item{pathway}{A single vector containing gene names.}

\item{groups}{A vector of length equal to the
number of rows in \code{x} containing two unique elements
(the two group names).}

\item{network_inference}{A function used to infer the pathway network.}

\item{n_perm}{The number of random permutations to perform during
permutation testing. If \code{n_perm == 1}, the permutation tests are not performed.}

\item{lp}{The lp value used to compute differential connectivity
scores. (Note: If a vector is provided, then the results are returned as
a list of \code{dnapath_list} objects, one result for each value of \code{lp}. This
option is available so that network inference methods only need to be run
once for each pathway when multple values of \code{lp} are being considered. This
may be useful when conducting simulation studies).}

\item{seed}{(Optional) Used to set.seed prior to permutation test for
each pathway. This allows results for individual pathways to be easily
reproduced.}

\item{verbose}{Set to TRUE to turn on messages.}

\item{...}{Additional arguments are passed into the network inference function.}
}
\value{
A list containing the results of the differential network analysis
for a single pathway.
}
\description{
This is an internal function used by \code{\link{dnapath}} to conduct
the differential network analysis on a single pathway.
}
\keyword{internal}
