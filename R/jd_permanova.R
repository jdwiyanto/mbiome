#' Calculates PERMANOVA
#'
#' @param ps A phyloseq object -- raw reads
#' @param var A variable(s) to test for PERMANOVA
#' @param strata A variable to stratify the PERMANOVA calculation. Defaults to NA
#' @return A list of PERMANOVA output
#' @export
jd_permanova = function(ps, var, strata = NA, mtd = 'euclidean') {

  sampledata = phyloseq::sample_data(ps)
  sampledata = subset(sampledata, !is.na(sampledata[[var]]))
  ps = phyloseq::phyloseq(sampledata, otu_table(ps), tax_table(ps))

  otu.perm <- phyloseq::otu_table(ps) %>% data.frame
  meta.perm <- phyloseq::sample_data(ps) %>% data.frame

  permanova <- vegan::adonis(t(otu.perm) ~ var %>% get,
                             strata = meta.perm[[strata]],
                             data = meta.perm,
                             permutations = 999,
                             method = mtd)

  return(permanova)
}
