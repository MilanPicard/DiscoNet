% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extraction_functions.R
\name{extract_subgraph_features}
\alias{extract_subgraph_features}
\title{Extract signature-based features}
\usage{
extract_subgraph_features(
  Graph,
  start_nodes = vnames(Graph),
  Signature_list,
  verbose = TRUE
)
}
\arguments{
\item{Graph}{A igraph object.}

\item{start_nodes}{Features are extracted for each of these nodes}

\item{verbose}{Default: TRUE. Wether or not to print the internal calculations}

\item{topo_results}{The results returned by the function extract_topolo}

\item{Signature}{A node signature (not list)}

\item{order}{The order of neighbors to consider}
}
\value{
Return something
}
\description{
Extract signature-based features from a biological network based on a list of signatures nodes.
}
