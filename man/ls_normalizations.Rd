% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.R
\name{ls_normalizations}
\alias{ls_normalizations}
\title{Library Size Normalization for SpatialExperiment Objects}
\usage{
ls_normalizations(SEO)
}
\arguments{
\item{SEO}{A `SpatialExperiment` object containing count data.}
}
\value{
The same `SpatialExperiment` object, with normalized counts added to
the 'logcounts' slot.
}
\description{
Normalizes the count data in a `SpatialExperiment` object based on library size.
This function calculates the library sizes (total counts per spot), computes the
mean library size, and scales the counts by these size factors. The normalized
counts are then stored in the 'logcounts' slot of the SpatialExperiment object.
}
\examples{
# Assuming 'SEO' is a SpatialExperiment object
SEO_normalized = ls_normalizations(SEO)

}
