#' Plots alpha-diversity with statistics
#'
#' @param ps A phyloseq object
#' @param var A variable(s) to classify visualisation
#' @param measure Alpha-diversity measures to visualise. Singular value. Defaults to Shannon
#' @param stat Statistical method to conduct. Defaults to Kruskal-Wallis
#' @export
jd_plot_alpha2 = function (ps, var, measure = "Shannon", stat = "kruskal.test")
{
  figlist = lapply(var, function(var) {
    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, !is.na(sampledata[[var]]))
    ps = phyloseq::phyloseq(otu_table(ps), sampledata, tax_table(ps))
    rich = phyloseq::estimate_richness(ps)
    rich$Pielou = rich$Shannon / log(rich$Observed)
    sampledata = cbind(sampledata, rich)
    fig1 = ggplot2::ggplot(sampledata, aes_string(x = var, y = measure,
                                                  fill = var)) + ggplot2::geom_boxplot() +
      ggpubr::stat_compare_means(method = stat) +
      ggplot2::theme_bw() + ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'blank') +
      ggplot2::ylab(measure) +
      ggplot2::xlab(stringr::str_to_sentence(var)) +
      ggplot2::geom_jitter(width = 0.1)
    return(fig1)
  })
  if (length(var) > 1) {
    fig = ggpubr::ggarrange(plotlist = figlist)
    return(fig)
  }
  else {
    fig = ggpubr::ggarrange(plotlist = figlist)
    return(fig)
  }
}

