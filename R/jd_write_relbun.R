#' Calculates the relative abundance value of taxa based on the chosen variable
#'
#' @param ps A phyloseq object -- raw reads
#' @param var Variable(s) to calculate relative abundance data
#' @return CSV files containing the relative abundance of each taxa per variable
#' @export
jd_write_relbun = function(ps, var) {

  for (ps in ps) {

    edge.filt <- phyloseq::otu_table(ps %>% get)
    threshold <- ( edge.filt %>% colSums %>% max ) * 0.01                           # filter taxa > 1% abundance
    keep <- ( edge.filt %>% colSums ) > threshold
    edge.filt2 <- edge.filt[,keep]
    ps.filt <- phyloseq::phyloseq(edge.filt2, sample_data(ps %>% get), tax_table(ps %>% get))


    rel.bun <- ps.filt %>%
      phyloseq::transform_sample_counts(function(x) {x/sum(x)} ) %>%
      phyloseq::psmelt()

    test <- subset(rel.bun, !is.na(rel.bun[['Abundance']]))

    for (var2 in var) {

      for (y in test$OTU %>% unique) {

        test2 <- subset(test, test$OTU == y)
        dftest <- tapply(test2$Abundance, test2[[var2]], mean) %>% as.data.frame
        colnames(dftest) <- y

        assign(paste0('rel.bun.', ps, '_', var2, '_', y), dftest)

      }

      testlist <- ls(pattern = paste0('rel.bun.', ps, '_', var2))
      testcomp <- dplyr::bind_cols(mget(testlist))

      assign(paste0('df.rel_', ps, '_', var2), testcomp)
      write.csv(assign(paste0('df.rel_', ps, '_', var2), testcomp), paste0('df.rel_', ps, '_', var2, '.csv'))
      write.csv(sample_data(ps %>% get), paste0('df.rel_', ps, '.meta.csv'))

    }
  }

}
