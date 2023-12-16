#' Library Size Normalization for SpatialExperiment Objects
#'
#' Normalizes the count data in a `SpatialExperiment` object based on library size.
#' This function calculates the library sizes (total counts per spot), computes the
#' mean library size, and scales the counts by these size factors. The normalized
#' counts are then stored in the 'logcounts' slot of the SpatialExperiment object.
#'
#' @param SEO A `SpatialExperiment` object containing count data.
#'
#' @return The same `SpatialExperiment` object, with normalized counts added to
#' the 'logcounts' slot.
#'
#' @examples
#' # Assuming 'SEO' is a SpatialExperiment object
#' SEO_normalized = ls_normalizations(SEO)
#'
#' @export
#' @import SingleCellExperiment
#' @import MatrixGenerics
ls_normalizations <- function(SEO) {
  if (!inherits(SEO, 'SpatialExperiment')) {
    stop("Input is not a SpatialExperiment object.")
  }

  # Calculate library sizes (sum of counts per spot)
  ls <- colSums(counts(SEO))

  # Calculate mean library size
  mean_ls <- mean(ls)

  # Calculate size factors
  size_fac <- ls / mean_ls

  # Apply normalization
  normalized_counts <- log1p(sweep(counts(SEO), 2, size_fac, FUN = "/"))

  # Update the counts in the SpatialExperiment object's 'logcounts' slot
  assay(SEO, "logcounts") <- normalized_counts

  return(SEO)
}
