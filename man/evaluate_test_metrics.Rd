% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneSelectR.R
\name{evaluate_test_metrics}
\alias{evaluate_test_metrics}
\title{Evaluate Test Metrics for a Grid Search Model}
\usage{
evaluate_test_metrics(grid_search, X_test, y_test, modules)
}
\arguments{
\item{grid_search}{A grid search object containing the best estimator.}

\item{X_test}{A data frame or matrix of test features.}

\item{y_test}{A vector of test labels.}

\item{modules}{A list of Python modules used in the function.}
}
\value{
A list containing key performance metrics of the best model:
- @field precision: The weighted precision score.
- @field recall: The weighted recall score.
- @field f1: The weighted F1 score.
- @field accuracy: The overall accuracy score.
These metrics are crucial for evaluating the effectiveness of the model on test data.
}
\description{
This function takes a grid search object, test data, and test labels to evaluate the performance
of the best model found during grid search.
}
