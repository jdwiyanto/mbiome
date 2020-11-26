#' Calculates alpha-diversity statistics
#'
#' @param ps A phyloseq object
#' @param var A variable to associate with the chosen alpha-diversity index
#' @param measure An alpha-diversity index
#' @export
jd_alpha_stat = function(ps, var, measure) {

  list = lapply(var, function(var) {

    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, !is.na(sampledata[[var]]))
    ps = phyloseq::phyloseq(otu_table(ps), sampledata, tax_table(ps))

    rich = phyloseq::estimate_richness(ps, measures = measure)
    data = sample_data(ps) %>% data.frame
    data[[measure]] = rich[[measure]]

    formula = as.formula(paste0(measure, ' ~ ', var))
    statistic = ggpubr::compare_means(formula, data = data)
    return(statistic)

  })

  return(list)

}
