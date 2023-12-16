#' Comprehensive Analysis Pipeline for SpatialExperiment Data
#'
#' This function provides an end-to-end analysis pipeline for `SpatialExperiment` objects.
#' It first normalizes the data using library size normalization, then selects the top
#' variable genes, and finally fits Negative Binomial models to the gene expression data.
#' This pipeline is designed to be a convenient wrapper for a series of complex steps,
#' streamlining the analysis of spatial transcriptomics data.
#'
#' @param seo A `SpatialExperiment` object to be processed.
#' @param n The number of top variable genes to select for modeling.
#'
#' @return A list of `glm.nb` model objects for each gene where the model fitting
#' was successful.
#'
#' @examples
#' # Assuming 'seo' is a preprocessed SpatialExperiment object
#' models <- modelFit(seo, n = 500)
#'
#' @export
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom scran modelGeneVar getTopHVGs
#' @importFrom MASS glm.nb
modelFit=function(seo,n){
  data=ls_normalizations(seo)
  data=select_features(data, n)
  fitNB(data)
}
