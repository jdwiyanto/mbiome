#' Calculates differential abundance statistics based on DESeq2
#'
#' @param ps A phyloseq object
#' @param var A variable to analyse
#' @param threshold A number to filter features with less than the chosen number log-fold changes. Defaults to 1.5
#' @export
jd_deseq = function(ps, var, threshold = 1.5) {

  loop = lapply(var, function(var) {
    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, !is.na(sampledata[[var]]))
    ps = phyloseq::phyloseq(otu_table(ps) + 1, sampledata, tax_table(ps))

    formula = as.formula(paste0('~ ', var))
    dds = phyloseq::phyloseq_to_deseq2(ps, formula) %>% DESeq2::DESeq()

    trial = lapply(2:length(DESeq2::resultsNames(dds)), function(x) {

      test = DESeq2::results(dds, name = DESeq2::resultsNames(dds)[x], alpha = 0.05, tidy = T) %>%
        subset(padj < 0.05) %>% data.frame(row.names = 1) %>%
        subset(abs(log2FoldChange) > threshold)  %>% tibble::rownames_to_column('taxa') %>% dplyr::arrange(-abs(log2FoldChange))

      if (test %>% nrow < 1) {

      } else {
        test$comparison = DESeq2::resultsNames(dds)[x]

        return(test)

      }
    })
    trialist = list(trial)
    names(trialist) = var
    return(trialist)
  })
  return(loop)
}
