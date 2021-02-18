#' Plots the relative abundance of a microbiome dataset
#'
#' @param ps A phyloseq object -- raw read counts
#' @param var A character string denoting the variable to be plotted. supports multiple variables
#' @param taxrank (for jd_plot_rel) A taxa level to visualise
#' @return A ggplot2 object
#' @export
jd_plot_rel <- function(ps, var, taxrank) {

  colnames(tax_table(ps)) = tolower(colnames(tax_table(ps)))

  figlist = lapply(var, function(xx) {

    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, !is.na(sampledata[[xx]]))
    ps2 = phyloseq::phyloseq(sampledata, otu_table(ps), tax_table(ps))

    if (taxrank == 'phylum') {

      p.rel <- ps2 %>% phyloseq::tax_glom(taxrank) %>%
        phyloseq::transform_sample_counts(function(x) {x/sum(x)} ) %>%
        phyloseq::psmelt() %>%
        dplyr::arrange(taxrank %>% get)

    } else {

      keep <-  (ps2 %>% phyloseq::taxa_sums %>% sort(decreasing = T))[1:20] %>% names
      p.rel <- phyloseq::prune_taxa(keep, ps2)
      p.rel <- p.rel %>% phyloseq::transform_sample_counts(function(x) {x/sum(x)} ) %>%
        phyloseq::psmelt() %>%
        dplyr::arrange(taxrank %>% get)

    }


    figrel = ggplot2::ggplot(p.rel, aes_string(x = xx, y = 'Abundance', fill = taxrank)) +
      ggplot2::geom_bar(stat = "identity", position = "fill") +
      ggplot2::labs(x = 'group', y = 'Relative abundance (%)') +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank(), legend.position = 'bottom') +
      ggplot2::guides(fill = guide_legend(nrow = 5, byrow = T)) +
      ggsci::scale_fill_d3(palette = 'category20')

    return(figrel)

  })

  if (var %>% length > 1) {

    fig = ggpubr::ggarrange(plotlist = figlist, common.legend = T, labels = 'auto')
    return(fig)

  } else {

    fig = ggpubr::ggarrange(plotlist = figlist, common.legend = T)
    return(fig)
  }
}

#' @export
jd_biome_rel = function(ps = ps.core.genus, testvar = testvar) {

  p.rel <- ps.genus.core %>% phyloseq::tax_glom('Phylum') %>%
    phyloseq::transform_sample_counts(function(x) {x/sum(x)} ) %>%
    phyloseq::psmelt() %>%
    dplyr::arrange(Phylum)

  test = p.rel[,c(testvar, 'Phylum', 'Abundance')]
  #test = na.omit(test)

  testdata = lapply(testvar, function(x) {
    data = test[, c(x, 'Phylum', 'Abundance')]
    data[[x]] = paste0(x, '_', tolower(data[[x]]))
    names(data) = c('var', 'Phylum', 'Abundance')
    return(data)
  })
  names(testdata) = testvar

  testdata2 = bind_rows(testdata)
  testdata2$Phylum = factor(testdata2$Phylum, levels = c('Proteobacteria', 'Desulfobacterota', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes'))

  fig = ggplot(testdata2, aes(x = var, fill = Phylum, y = Abundance)) + geom_bar(stat = 'identity', position = 'fill') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab(element_blank())
  return(fig)
}
