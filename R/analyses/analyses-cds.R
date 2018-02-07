analyses$CDS_validation <- function(given, save_as, genome) {
  
  n_stop_codons <- function(x) {
    str_locate_all(x, 'TGA|TAA|TAG') %>% map_int(~sum(.[, 'start'] %% 3 == 1L))
  }
  
  cds <- read_csv(given, col_types = cols())
  
  valid <- 
    cds %>%
    mutate(exon_seq = get_genomic_sequence(chr, strand, start, end, genome)) %>%
    group_by(tx) %>%
    arrange(tx, exon) %>%
    summarise(dna = str_c(exon_seq, collapse = '')) %>%
    ungroup() %>%
    filter(
      near(0, str_length(dna) %% 3),     # length is mutiple of 3
      str_detect(dna, '^ATG'),           # begins with start
      near(1, n_stop_codons(dna)),       # only one stop codon
      str_detect(dna, 'TGA$|TAA$|TAG$')  # stop codon is at end
    )
  
  filter(cds, tx %in% valid$tx) %>%
    write_csv(save_as)
}

analyses$CDS_validation_hg38 <- function(given   = 'data/CDS/Hsapiens-UCSC-hg38.csv',
                                         save_as = 'data/CDS/Hsapiens-UCSC-hg38-validated.csv',
                                         genome  = BSgenome.Hsapiens.UCSC.hg38::Hsapiens) {
  analyses$CDS_validation(given, save_as, genome)
}

analyses$CDS_validation_hg19 <- function(given   = 'data/CDS/Hsapiens-UCSC-hg19.csv',
                                         save_as = 'data/CDS/Hsapiens-UCSC-hg19-validated.csv',
                                         genome  = BSgenome.Hsapiens.UCSC.hg19::Hsapiens) {
  analyses$CDS_validation(given, save_as, genome)
}
