% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extraction_functions.R
\name{DEG_in_neighborhood}
\alias{DEG_in_neighborhood}
\title{Extract signature-based features}
\usage{
DEG_in_neighborhood(Graph, start_nodes, degs_list, order = 1, verbose = TRUE)
}
\arguments{
\item{Graph}{A igraph object.}

\item{start_nodes}{Features are extracted for each of these nodes}

\item{order}{The order of neighbors to consider}

\item{degs}{A node signature (not list)}
}
\value{
Return something
}
\description{
Extract signature-based features from a biological network based on a list of signatures nodes.
}
