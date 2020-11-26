#' Spiec Easi network analysis
#'
#' @param ps A Phyloseq object
#' @param var A variable(s) to test for
#' @return list of Spiec Easi object and ggplot2 object
#' @export
jd_spiec = function(ps, var) {

  all = lapply(var, function(var) {

    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, country == var)
    ps2 = phyloseq::phyloseq(sampledata, otu_table(ps), tax_table(ps))

    se<- SpiecEasi::spiec.easi(ps2, method='mb', lambda.min.ratio=1e-2,
                               nlambda=20, pulsar.params=list(rep.num=50))

    ig2 <- SpiecEasi::adj2igraph(getRefit(se),  vertex.attr=list(name=taxa_names(ps2)))

    fig = phyloseq::plot_network(ig2, ps2, type='taxa', color="phylum") + theme(legend.position = 'bottom')

    output = list(se, ig2, fig)

    return(output)

  })

  return(all)

}
