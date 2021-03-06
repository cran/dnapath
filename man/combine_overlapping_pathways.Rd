% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reactome.R
\name{combine_overlapping_pathways}
\alias{combine_overlapping_pathways}
\title{Modify a pathway list to combine overlapping pathways.}
\usage{
combine_overlapping_pathways(pathway_list, overlap_limit = 0.9)
}
\arguments{
\item{pathway_list}{A list of pathways obtained from
\code{\link{get_reactome_pathways}}.}

\item{overlap_limit}{A percentage between 0 and 1. If two pathways
overlap by more than this amount, they are combined into one pathway.}
}
\value{
A modified list with overlapping pathways combined together.
}
\description{
Modify a pathway list to combine overlapping pathways.
}
\keyword{internal}
