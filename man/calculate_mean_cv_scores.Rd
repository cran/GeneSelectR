% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneSelectR.R
\name{calculate_mean_cv_scores}
\alias{calculate_mean_cv_scores}
\title{Calculate Mean Cross-Validation Scores for Various Feature Selection Methods}
\usage{
calculate_mean_cv_scores(selected_pipelines, cv_best_score)
}
\arguments{
\item{selected_pipelines}{A list of pipelines for different feature selection methods.}

\item{cv_best_score}{A list or vector of cross-validation scores.}
}
\value{
A dataframe containing the mean and standard deviation of cross-validation scores for each method.
}
\description{
Calculate Mean Cross-Validation Scores for Various Feature Selection Methods
}
