#' Plot alpha diversity based on selected variables
#'
#' @param ps A phyloseq object
#' @param var Variable(s) to organise visualisation
#' @param measures Alpha-diversity measures (for jd_plot_alpha, defaults to Chao1 and Shannon. For jd_plot_alpha_stat, defaults to Shannon. Singular value only)
#' @param stat Statistical method to conduct. Defaults to Kruskal-Wallis
#' @export
jd_plot_alpha <- function(ps, var, measures = c('Chao1', 'Shannon')) {


  figlist = lapply(var, function(var) {
    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, !is.na(sampledata[[var]]))
    ps2 = phyloseq::phyloseq(sampledata, otu_table(ps), tax_table(ps))

    fig = phyloseq::plot_richness(ps2,
                                  x = var,
                                  measures = measures) +
      ggplot2::geom_violin(fill = 'white', alpha = 0.5) +
      ggplot2::geom_boxplot(alpha = 1) +
      ggplot2::theme_bw() +
      ggplot2::xlab(element_blank()) +
      ggplot2::scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))

    return(fig)

  })

  if( var %>% length > 1 ) {

    fig = gggpubr::ggarrange(plotlist = figlist, labels = 'auto')

    return(fig)
  } else {

    fig = ggpubr::ggarrange(plotlist = figlist)

    return(fig)
  }

}

#' @export
jd_plot_alpha_stat = function (ps, var, measure = "Shannon", stat = "kruskal.test")
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

