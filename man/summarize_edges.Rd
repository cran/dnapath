% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summarize_edges}
\alias{summarize_edges}
\title{Summarize differential connections for a pathway}
\usage{
summarize_edges(x, alpha = 0.1, monotonized = FALSE, require_dc_genes = FALSE)
}
\arguments{
\item{x}{A 'dnapath' object from \code{\link{dnapath}}.}

\item{alpha}{Threshold for p-values of edge DC scores.
Defaults to 0.1 or the minimum possible threshold for the number
of permutations performed, whichever is greater.}

\item{monotonized}{If TRUE, monotonized p-values are used.}

\item{require_dc_genes}{If TRUE, the gene-level differential connectivity
p-value of the two genes for a given edge are also considered when deciding
whether an edge is differentially connected. If neither gene is significantly
differentially connected, then the edge between them will not be either.}
}
\value{
A tibble summarizing the differential connections in
the pathway.
}
\description{
Summarize differential connections for a pathway
}
\examples{
data(meso)
data(p53_pathways)
set.seed(0)
results <- dnapath(x = meso$gene_expression, pathway_list = p53_pathways,
                   group_labels = meso$groups, n_perm = 10)
summarize_edges(results[[1]])
}
\seealso{
\code{\link{summarize_pathways}}, \code{\link{summarize_genes}}
}
