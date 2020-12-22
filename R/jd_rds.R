#' Merge ASV sequence tables (RDS files) from DADA2 and annotate taxonomy using SILVA
#'
#' Function must be run on main folder containing all seqtab.RDS files to be annotated (recursive folders supported).
#' This function will output abundance and taxonomy tables in the working directory.
#'
#' @param seqlength_min minimum sequence length to keep, default to 400 for V3-V4 region
#' @param seqlength_max maximum sequence length to keep, default to 470 for V3-V4 region
#' @param refpath path to SILVA reference database
#' @return csv files of abundance and taxonomy table
#' @export
jd_rds = function(seqlength_min = 400, seqlength_max = 470,
                  refpath = "E:/ubuntu/virtual_shared/silva_nr99_v138_wSpecies_train_set.fa.gz") {

  tbm <- dir(pattern = '.rds$', recursive = T)
  stm <- dada2::mergeSequenceTables(tables = tbm)
  saveRDS(stm, 'seqtab_compiled.rds')

  #remove chimera
  stm.nochim <- dada2::removeBimeraDenovo(stm, method="consensus", multithread=TRUE, verbose=TRUE)
  saveRDS(stm.nochim, 'seqtab_compiled_nochim.rds')

  # filter stm to remove unwanted sequence length
  stm.nochim2 <- subset(stm.nochim,
                        select = dplyr::between(nchar(dada2::getSequences(stm.nochim)), seqlength_min, seqlength_max))
  saveRDS(stm.nochim2, 'seqtab_compiled_nochim_v34_filtered.rds')

  # assign taxonomy
  taxa <- dada2::assignTaxonomy(stm.nochim2, refpath, multithread=TRUE)

  taxa.print <- taxa # Removing sequence rownames for display only
  rownames(taxa.print) <- NULL

  ps <- phyloseq::phyloseq(phyloseq::otu_table(stm.nochim2, taxa_are_rows=FALSE), phyloseq::tax_table(taxa))

  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- phyloseq::merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)), "_", tax_table(ps)[,6])

  write.csv(phyloseq::tax_table(ps), 'phyloseq_tax.csv')
  write.csv(phyloseq::otu_table(ps), 'phyloseq_otu.csv')

}

