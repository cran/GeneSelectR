#' @title Class Union for Test metrics output that could contain either a dataframe or a lists
#' @description A class union that can contain either a data frame or a list.
#' @docType class
#' @seealso \link{setClassUnion}
#' @keywords internal
#' @export
setClassUnion("TestMetrics", c("data.frame", "list"))

#' PipelineResults class
#'
#' A class to hold the results of GeneSelectR function.
#'
#' @slot best_pipeline A named list containing parameters of the best performer pipeline
#' @slot cv_results A list of the cross-validation results for each pipeline.
#' @slot inbuilt_feature_importance A list of the inbuilt mean feature importances for each method across all splits.
#' @slot permutation_importance A list of the permutation importances for each method.
#' @slot cv_mean_score A data.frame of mean scores from cross-validation.
#' @slot test_metrics A data.frame containing metrics (F1, accuracy, precision, and recall) calculated on the unseen test set. Contains mean values across splits as well as standard deviation.
#' @importFrom methods representation
#' @exportClass PipelineResults
setClass("PipelineResults",
         slots = list(
           best_pipeline = "list",
           cv_results = "list",
           inbuilt_feature_importance = "list",
           permutation_importance = 'list',
           cv_mean_score = 'data.frame',
           test_metrics = 'TestMetrics'
         ))

#' @title GeneList class
#' @description A class to hold annotated gene list for a single method.
#' @slot SYMBOL A character vector of gene names.
#' @slot ENSEMBL A character vector of Ensembl IDs.
#' @slot ENTREZID A character vector of Entrez IDs.
#' @exportClass GeneList
setClass("GeneList",
         slots = list(
           SYMBOL = "character",
           ENSEMBL = "character",
           ENTREZID = "character"
         ))


#' @title AnnotatedGeneLists class
#' @description A class to hold a list of GeneList objects, each representing a method.
#' @slot inbuilt A list of GeneList objects containing annotations for genes selected with inbuilt, model-specific feature importance.
#' @slot permutation A list of GeneList objects containing annotations for genes selected with permutation importance.
#' @exportClass AnnotatedGeneLists
setClass("AnnotatedGeneLists",
         methods::representation(
            inbuilt = "list",
            permutation = 'list'
         ))

