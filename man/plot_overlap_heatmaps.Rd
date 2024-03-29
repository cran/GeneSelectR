% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_overlap_heatmaps.R
\name{plot_overlap_heatmaps}
\alias{plot_overlap_heatmaps}
\title{Generate Heatmaps to Visualize Overlap and Similarity Coefficients between Feature Lists}
\usage{
plot_overlap_heatmaps(coefficients, save_plot = FALSE, filename = NULL)
}
\arguments{
\item{coefficients}{A list of matrices showing the Overlap, Jaccard, and Soerensen-Dice coefficients for the feature lists.}

\item{save_plot}{A logical value indicating whether to save the heatmap plots to a file or not. Default is FALSE.}

\item{filename}{A character string specifying the filename for the saved heatmap plots (if save_plot = TRUE).}
}
\value{
A grid of ggplot2 heatmap objects visualizing the Overlap, Jaccard, and Soerensen-Dice coefficients.
The grid layout includes heatmaps for both 'inbuilt' and 'permutation' feature importance coefficients (if available).
If \code{save_plot} is TRUE, the heatmaps are also saved to the specified file.
}
\description{
This function takes a list of matrices of overlap and similarity coefficients and generates heatmaps to visualize them.
}
\examples{
\donttest{
# Assuming `coefficients` is a list containing matrices for Overlap, Jaccard,
# and Soerensen-Dice coefficients
# For demonstration, let's create a mock coefficients list
mock_matrix <- matrix(runif(25), nrow = 5)
coefficients <- list(inbuilt_feature_importance_coefficient = list(overlap = mock_matrix,
                     jaccard = mock_matrix, soerensen = mock_matrix),
                     permutation_importance_coefficients = list(overlap = mock_matrix,
                     jaccard = mock_matrix, soerensen = mock_matrix))

# Plot the overlap heatmaps
heatmap_plots <- plot_overlap_heatmaps(coefficients)
print(heatmap_plots)
}
}
