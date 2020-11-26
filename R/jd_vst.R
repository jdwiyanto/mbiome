#' Transform raw read counts using VariableStabilizingTransformation in DESeq2
#'
#' @param ps A phyloseq object
#' @return A phyloseq object with a VST-transformed reads
#' @export
jd_vst = function(ps) {

  pseudo = phyloseq::otu_table(ps) + 1
  ps2 = phyloseq::phyloseq(pseudo, tax_table(ps), sample_data(ps))
  dds = phyloseq::phyloseq_to_deseq2(ps2, ~1)

  dds.vst <- DESeq2::varianceStabilizingTransformation(dds, blind = F) # transform with VST
  dds.vst <-  assay(dds.vst)
  dds.vst[dds.vst < 0.0] <-  0.0 # convert negative value to 0 (for distance matrix analysis)

  ps.vst = phyloseq::phyloseq(otu_table(dds.vst, taxa_are_rows = T), sample_data(ps), tax_table(ps))

  return(ps.vst)

}
