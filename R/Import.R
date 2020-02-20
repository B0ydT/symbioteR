#' Read QIIME 2 Output from a Folder
#'
#' @param filepath A filepath
#'
#' @return A phyloseq object
#' @export
#' @importFrom qiime2R qza_to_phyloseq read_qza
#' @importFrom magrittr %>% set_colnames set_rownames extract2
#' @importFrom stringr str_replace_all
#' @importFrom utils getFromNamespace
#' @import qpcR
#'
readQIIMEFolder <- function(filepath) {
  rbind.na = getFromNamespace("rbind.na", "qpcR")
  physeq <- qza_to_phyloseq(
    features = paste0(filepath, "/feat_tab.qza"),
    #tree = paste0(filepath, "/phylogeny.qza"),
    taxonomy =  paste0(filepath, "/taxonomy.qza"),
    metadata = paste0(filepath, "/metadata.txt"))
  tax <- paste0(filepath, "/taxonomy.qza") %>%
    read_qza()
  tax_table(physeq) <- as.character(tax$data$Taxon) %>%
    strsplit(";") %>%
    do.call(what = rbind.na) %>%
    set_colnames(c("Kingdom", "Phylum", "Class", "Order", "Family",
                   "Genus", "Species")) %>%
    data.frame() %>%
    lapply(str_replace_all, "^._._*", "") %>%
    as.data.frame() %>%
    set_rownames(tax$data$Feature.ID) %>%
    as.matrix()
  physeq <- paste0(filepath, "/rep_seqs.qza") %>%
    read_qza() %>%
    extract2("data") %>%
    merge_phyloseq(physeq)
  return(physeq)
}
