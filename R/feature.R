#' Select Top Variable Genes from a SpatialExperiment Object
#'
#' @param SEO A `SpatialExperiment` object to analyze.
#' @param n The number of top variable genes to select.
#' @param assay_name The name of the assay to use for variability calculations.
#'
#' @return A subset of the `SpatialExperiment` object containing only the top `n`
#' most variable genes.
#'
#' @examples
#' # Assuming 'SEO' is a SpatialExperiment object with logcounts
#' SEO_hvg = select_features(SEO, n = 500, assay_name = "logcounts")
#'
#' @export
#' @import SpatialExperiment
#' @import SummarizedExperiment
select_features <- function(SEO, n, assay_name = "logcounts") {
  if (!inherits(SEO, 'SpatialExperiment')) {
    stop("The provided object is not a Spatial Experiment object.")
  }

  if (!assay_name %in% assayNames(SEO)) {
    available_assays <- paste(assayNames(SEO), collapse = ", ")
    stop("Specified assay not found in the SpatialExperiment object. Available assays are: ", available_assays, ".")
  }

  # Compute model var. of gene's expression across spatial locations using the specified assay
  gene_var <- scran::modelGeneVar(SEO, assay.type = assay_name)

  # Select top n most variable genes
  hvg <- scran::getTopHVGs(gene_var, n = n)

  # Subset the SpatialExperiment object
  spe_hvg <- SEO[hvg, ]

  return(spe_hvg)
}
