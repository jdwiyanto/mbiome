#' Transform raw read counts using centered-log transformation approach
#'
#' @param ps a phyloseq object with raw read counts. taxa have to be the columns (ensure that otu_table(ps, taxa_are_rows = F))
#' @return a phyloseq object with CLR-transformed read counts (with taxa_are_rows = T)
#' @export
jd_clr <- function(ps) {

  edge.propr <- phyloseq::otu_table(ps)
  edge.propr2 <- edge.propr %>% as.data.frame() %>% t() %>% propr::propr()
  edge.propr3 <- edge.propr2@logratio %>% phyloseq::otu_table(taxa_are_rows = T)

  ps.clr <- phyloseq::phyloseq(edge.propr3, sample_data(ps), tax_table(ps))

  return(ps.clr)

}
