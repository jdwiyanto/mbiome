#' Extract taxa names
#'
#' @param ps A phyloseq object from where the taxa names will be extracted
#' @param df A dataframe, whose rownames are taxa codes
#' @return data frame containing taxonomy names
#' @export
jd_annotate = function(ps, df) {
  toget = df %>% rownames
  taxname = phyloseq::tax_table(ps)[toget,] %>% data.frame
  ordname = taxname$order
  famname = taxname$family
  genname = taxname$genus
  spname = taxname$species
  txname = taxname$taxon
  data = data.frame(order = ordname, family = famname, genus = genname, species = spname, taxon = txname)
  return(data)
}

# return taxa name for pathway ps

jd_annotate_path = function(ps, df) {
  toget = df %>% rownames
  taxname = phyloseq::tax_table(ps)[toget,] %>% data.frame
  pathname = taxname$pathway
  data = data.frame(pathway = pathname)
  return(data)
}
