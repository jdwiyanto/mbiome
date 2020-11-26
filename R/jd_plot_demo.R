#' Plot demographic distribution
#'
#' This functions plots the demographic distribution of the chosen variable(s)
#'
#' @param ps A phyloseq object
#' @param var A variable to plot
#' @param fill A variable to classify
#' @return a ggplot2 object
#' @export
jd_plot_demo <- function(ps, var, fill) {

sampledata <- sample_data(ps)

figlist = lapply(var, function(var) {

  figure <- ggplot2::ggplot(sampledata[!is.na(sampledata[[var]]),], aes(x = var %>% get, fill = fill %>% get)) +
    ggplot2::geom_bar(stat = 'count') +
    ggplot2::xlab(var) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())

  assign(paste0('fig.demo.', var), figure)
})

if( var %>% length > 1 ) {

  fig = ggpubr::ggarrange(plotlist = figlist, labels = 'auto')

  return(fig)

} else {

  fig = ggpubr::ggarrange(plotlist = figlist)

  return(fig)
}
}

