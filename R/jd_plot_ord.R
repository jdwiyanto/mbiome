#' Ordinates metagenomic data based on variable(s)
#'
#' @param ps A phyloseq object
#' @param var A variable to ordinate
#' @param method Ordination plot used
#' @param distance A distance matrix for the ordination
#' @export
jd_plot_ord <- function(ps, var, method = 'PCoA', distance ='euclidean') {

figlist = lapply(var, function(var) {
  sampledata = phyloseq::sample_data(ps)
  sampledata = subset(sampledata, !is.na(sampledata[[var]]))
  ps2 = phyloseq::phyloseq(sampledata, otu_table(ps), tax_table(ps))

  ord <- phyloseq::ordinate(ps2, method = method, distance = distance)

  fig <- phyloseq::plot_ordination(ps2,
                                   type = "samples",
                                   ordination = ord,
                                   color = var) +
    ggplot2::theme_bw() +
    ggsci::scale_color_d3()

  return(fig)

})

if( var %>% length > 1 ) {

  fig = ggpubr::ggarrange(plotlist = figlist, labels = 'auto')

  return(fig)

} else {

  fig = ggpubr::ggarrange(plotlist = figlist)

  return(fig)
}
}
