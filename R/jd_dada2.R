#' Run DADA2 ASV inference from primer- and adapter-removed fastq files
#'
#' Function must be run on folders containing the fastq files to be inferred to ASV. This function will output several files,
#' including the list of filtered reads, sequence table file in RDS form (chimera not removed!), redundant fasta files for PAPRICA input,
#' and save the whole environment (.R file)
#'
#' It is important to take note that different sample batch should be run separately. Recommended to keep each batch on separate folder
#' and run this function on a folder loop. For example:
#'
#' for (directory in 1:length(list.dirs(recursive = F))) {
#' path = getwd()
#' dire = list.dirs(recursive = F)
#' setwd(dire[directory])
#' jd_dada2()
#' setwd(path)
#' }
#'
#'
#' @param seqlength_min minimum sequence length to keep, default to 400 for V3-V4 region
#' @param seqlength_max maximum sequence length to keep, default to 470 for V3-V4 region
#' @param trunclen = dada2 input determining the length of sequence to keep after quality control. Default to no trimming c(0,0)
#' @return csv files of filtered reads summary, ASV sequence table in RDS format, PAPRICA-friendly redundant fasta sequence files, saved R environment
#' @export
jd_dada2 = function(seqlength_min = 400, seqlength_max = 470, truncLen = c(0,0)) {

  path = getwd()

  fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))           #set pattern to file names
  fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

  filtFs <- file.path(path, "dada2_filt", paste0(sample.names, "_1.fastq"))       # remove N from the reads - DADA can read no N - nucleotide
  filtRs <- file.path(path, "dada2_filt", paste0(sample.names, "_2.fastq"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names

  out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = truncLen, maxEE = c(2,6), maxN = 0)
  write.csv(out, 'filtered_reads.csv')

  errF <- dada2::learnErrors(filtFs, multithread=TRUE, verbose = 1)                      # learn error rates
  errR <- dada2::learnErrors(filtRs, multithread=TRUE, verbose = 1)

  dadaFs <- dada2::dada(filtFs, err=errF, multithread = TRUE)                            # sample inferences
  dadaRs <- dada2::dada(filtRs, err=errR, multithread = TRUE)

  #merge read
  mergers <- dada2::mergePairs(dadaFs, filtFs, dadaRs, filtRs, maxMismatch = 0, verbose = TRUE)
  seqtab <- dada2::makeSequenceTable(mergers)

  dir.create(paste0(path, '/', 'merged/'))                                       # create redundant seq files for paprica

  for(name in names(mergers)){
    temp <- mergers[[name]]
    temp <- temp[which(nchar(temp$sequence) %in% seqlength_min:seqlength_max),]
    write.csv(temp, paste0(path, '/merged/', name, '.csv'), quote = F, row.names = F)
  }

  saveRDS(seqtab, 'seqtab.rds')                                                   # write to file
  save.image('dada2_environment.R')

}
