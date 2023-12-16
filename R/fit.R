#' Fit Negative Binomial Models to Gene Expression Data
#'
#' Fits Negative Binomial models to gene expression data contained in a
#' `SpatialExperiment` object. This function iterates over each gene, fits a model
#' using gene expression counts, and stores the resulting model if successful.
#'
#' @param data A `SpatialExperiment` object containing gene expression data.
#' @param assay_type The type of assay to use from the SpatialExperiment object,
#' with the default being 'counts'.
#'
#' @return A list of `glm.nb` model objects for each gene, where the model fitting
#' was successful.
#'
#' @examples
#' # Assuming `top_500` is a preprocessed SpatialExperiment object
#' nb_models <- fitNB(top_500)
#'
#' @export
#' @import SpatialExperiment
#' @import MatrixGenerics
#' @import MASS
#' @importFrom stats glm.control

fitNB <- function(data, assay_type = 'counts') {
  # Extract counts matrix
  counts_matrix <- assay(data, assay_type)

  # Calculate library sizes and log-transform them
  library_sizes <- colSums(counts(data))
  logLS <- log(library_sizes)

  # Extract Annotated Cluster information
  AnnotatedCluster <- colData(data)$AnnotatedCluster

  # Initialize a list to store models
  nb_models <- vector("list", nrow(counts_matrix))

  # Iterate over each gene
  for (i in 1:nrow(counts_matrix)) {
    # Extract gene expression counts for the current gene
    gene_counts <- counts_matrix[i, ]

    # Create a dataframe for the model
    df <- data.frame(gene_expression = gene_counts,
                     AnnotatedCluster = factor(AnnotatedCluster),
                     logLS = logLS)

    # Check for non-finite values
    if (any(!is.finite(df$gene_expression)) || any(!is.finite(df$logLS))) {
      next
    }

    # Fit the Negative Binomial GLM with error handling
    nb_models[[i]] <- tryCatch({
      MASS::glm.nb(gene_expression ~ AnnotatedCluster + logLS, data = df, control = glm.control(maxit = 1000))
    }, error = function(e) {
      NULL  # Return NULL in case of an error
    })
  }

  # Count successful models
  successful_models <- sum(sapply(nb_models, function(x) !is.null(x)))

  cat('Successful models (without spatial data):', successful_models, "\n")

}
