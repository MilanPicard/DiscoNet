% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FS_and_ML.R
\name{InformationGain_perfect}
\alias{InformationGain_perfect}
\title{Internal function: calculate the maximum information gain for a feature}
\usage{
InformationGain_perfect(class)
}
\arguments{
\item{class}{A boolean vector representing the class (TRUE/FALSE)}
}
\value{
A single value representing the maximum information gain.
}
\description{
#' Calculate the maximum information gain for a dataset to later normalize the values obtained so that they go to 1.
}
