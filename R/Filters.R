#' Remove Contaminant Sequences from Bacterial Datasets
#'
#' @param physeq A phyloseq object.
#'
#' @return A phyloseq object.
#' @import phyloseq
#' @importFrom magrittr %>%
#' @export
#'
removeNonbacterial <- function(physeq) {
  Kingdom <- Order <- Family <- NULL
  subset_taxa(physeq, Kingdom != "Unassigned") %>%
    #Remove Taxa Unassigned at Kingdom Level
    subset_taxa(Order!= "Chloroplast" | is.na(Order)) %>% #Remove Chloroplasts
    subset_taxa(Family!= "Mitochondria" | is.na(Family))} #Remove Mitochondria
#' Apply a prevalence filter
#'
#' Reatin taxa observed in greater or equal to number of samples*threshold.
#'
#' @param physeq A phyloseq object.
#' @param threshold Prevalence threshold.
#'
#' @return A phyloseq object.
#' @import phyloseq
#' @export
#'
prevalenceFilter <- function(physeq, threshold) {
  physeq2 <- makeObservationsBinary(physeq)
  physeq2 <- tryCatch(prune_taxa(taxa_sums(
    physeq2) >= (threshold * nsamples(physeq)),
    physeq), error = function(e){
      return(NULL)
    })
  return(physeq2)
}
