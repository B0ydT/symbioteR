#' Convert an OTU Count Table to Presence/Absence
#'
#' @param physeq A phyloseq object.
#'
#' @return A phyloseq object.
#' @export
#'
makeObservationsBinary <- function(physeq) {
  physeq <- transform_sample_counts(physeq, function(x) {x/x})
  otu_table(physeq) <- transform_sample_counts(physeq, function(x) {
    replace_na(x, 0)})
  return(physeq)
}
#' Convert an OTU Count Table to Relative Abundance
#'
#' @param physeq A phyloseq object.
#'
#' @return A phyloseq object.
#' @export
#'
makeObservationsRelative <- function(physeq) {
  physeq <- transform_sample_counts(physeq, function(x) {x/sum(x)})
  return(physeq)
}
