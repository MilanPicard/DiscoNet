% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FS_and_ML.R
\name{AdaLASSO}
\alias{AdaLASSO}
\title{Adaptive LASSO for feature selection}
\usage{
AdaLASSO(df, class, nfolds = 10, verbose = TRUE)
}
\arguments{
\item{df}{A numerical dataframe with features as columns and observations as rows}

\item{class}{A boolean vector representing the class (TRUE/FALSE)}

\item{nfolds}{The number of folds during cross validation for tuning lambda}

\item{verbose}{Default to TRUE.Printing internal calculations.}
}
\value{
A string vector returning the selected features.
}
\description{
Return the list of features chosen through Adaptive LASSO during cross validation
}
\examples{
if(requireNamespace("mlbench", quietly = TRUE)) {
data("BreastCancer", package = "mlbench")
}

# Remove rows with missing values
BreastCancer = BreastCancer[rowSums(is.na(BreastCancer)) == 0, ]

# Remove the ID column and the class column.
df = BreastCancer[, -c(1, 11)]

# Convert each column to numeric class
df = as.data.frame(apply(df, 2, as.numeric))

# Retrieve the class as a logical vector
class = BreastCancer$Class == "malignant"

# Calculate selected features
BreastCancer_AdaLASSO = AdaLASSO(df, class, nfolds = 10)
}
