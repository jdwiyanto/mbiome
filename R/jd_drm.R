#' Calculates and plots Dirichlet Multinomial Mixture model of the dataset
#'
#' @param ps A phyloseq object -- raw read counts
#' @param nvar Number of variables to optimise number of enterotypes
#' @param var A variable to classify the detected enterotypes
#' @return A list of the detected enterotypes composition, metadata on the participants' enterotype, and ggplot2 object
#' @export
jd_drm <- function(ps, nvar, var) {

  d.otu <- ( otu_table(ps) ) + 1                                             # set otu table and add pseudo-count of 1
  cnts <- log10(colSums(d.otu))                                              # optional, confirm density plot is normal (neg binomial distribution)
  fit <- mclapply(1:nvar, dmn, count = d.otu, verbose = T)
  lplc <- sapply(fit, laplace)
  best <- fit[[which.min(lplc)]]

  d.type <- DirichletMultinomial::mixture(best)                              # give cluster prediction for each sample
  colnames(d.type) <- c(1:which.min(lplc))                                   # give column names to cluster prediction
  cluster <- colnames(d.type)[apply(d.type, 1, which.max)]                   # assign new column with chosen best prediciton
  d.type2 <- cbind(d.type, cluster) %>% as.data.frame                        # combine best cluster column with cluster prediction data frame

  p0 <- fitted(fit[[1]], scale=TRUE)                                         # scale by theta
  p3 <- fitted(best, scale=TRUE)
  colnames(p3) <- paste("m", 1:which.min(lplc), sep="")

  diff <- rowSums(abs(p3 - as.vector(p0)))                                   # generate dataframe of enterotypes
  o <- order(diff, decreasing=TRUE)
  cdiff <- cumsum(diff[o]) / sum(diff)
  df.drm <- ( head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 30) ) %>% as.data.frame()
  df.drm$names <- rownames(df.drm)
  df.drm2 <- df.drm[1:nvar,c((1:which.min(lplc) +1), (which.min(lplc) + 4))] %>%
    tidyr::pivot_longer(!names, names_to = 'enterotype', values_to = 'prop')

  fig.drm <- ggplot2::ggplot(df.drm2, aes( x = enterotype,
                                           y = prop*100,
                                           fill = names)) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::ylab('Proportion (%) explained by top 5 genera') +
    ggplot2::theme(legend.title = element_blank(), legend.position = 'bottom')


  for (x in var) {

    drm.table <- table(d.type2$cluster, df[[x]]) %>% prop.table(margin = 1) %>% data.frame

    fig.drm.prev <- ggplot2::ggplot(drm.table, aes(x = Var1,
                                                   y = Freq*100,
                                                   fill = Var2)) +
      ggplot2::geom_col(position = 'stack') +
      ggplot2::scale_fill_discrete(name = x) +
      ggplot2::theme_bw() +
      ggplot2::ylab('Proportion (%)') +
      ggplot2::xlab('Enterotype')

    assign(paste0('fig.drm.prev.', x), fig.drm.prev)

  }

  fig.drm.list <- paste0('fig.drm.prev.', var)
  fig.drm2 <- ggpubr::ggarrange(plotlist = mget(fig.drm.list), labels = 'auto')

  drmoutput <- list(which.min(lplc), df.drm, df.drm2, fig.drm, fig.drm2)
  names(drmoutput) <- c('best no of component', 'df of DRM', 'prop table of DRM', 'figure - top 5 genera', 'figure - enterotype distribution')

  return(drmoutput)

}
