#' Plot alpha diversity based on selected variables
#'
#' @param ps A phyloseq object
#' @param var Variable(s) to organise visualisation
#' @param measures Alpha-diversity measures (default to Chao1 and Shannon)
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
