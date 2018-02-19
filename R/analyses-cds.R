# ---- CDS validation ----
analysis_cds_validation_hg38 <- function(given   = 'data/CDS/Hsapiens-UCSC-hg38.csv',
                                         save_as = 'data/CDS/Hsapiens-UCSC-hg38-validated.csv',
                                         genome  = BSgenome.Hsapiens.UCSC.hg38::Hsapiens) {
  read_csv(given, col_types = cols()) %>%
    .validate_cds(genome) %>%
    write_csv(save_as)
}

analysis_cds_validation_hg19 <- function(given   = 'data/CDS/Hsapiens-UCSC-hg19.csv',
                                         save_as = 'data/CDS/Hsapiens-UCSC-hg19-validated.csv',
                                         genome  = BSgenome.Hsapiens.UCSC.hg19::Hsapiens) {
  read_csv(given, col_types = cols()) %>%
    .validate_cds(genome) %>%
    write_csv(save_as)
}

# ---- Transcript features ----
analysis_cds_tx_features_hg38 <- function(given   = 'data/CDS/Hsapiens-UCSC-hg38-validated.csv',
                                          save_as = 'data/CDS/Hsapiens-UCSC-hg38-validated-tx-features.csv',
                                          genome  = BSgenome.Hsapiens.UCSC.hg38::Hsapiens) {
  read_csv(given, col_types = cols()) %>%
    # Extract and construct coding sequenes from CDS exon coordinates
    mutate(exon_seq = get_genomic_sequence(chr, strand, start, end, genome)) %>%
    arrange(gene, tx, exon) %>%
    group_by(gene, tx, chr) %>%
    summarise(
      sequence   = str_c(exon_seq, collapse = ''),
      cds_length = str_length(sequence),
      aa_length  = cds_length / 3,
      n_exons    = n(),
      n_g        = str_count(sequence, 'G'),
      n_c        = str_count(sequence, 'C'),
      n_a        = str_count(sequence, 'A'),
      n_t        = str_count(sequence, 'T')
    ) %>%
    select(-sequence) %>%
    group_by(gene, chr) %>%
    mutate(n_isoforms = n()) %>% # isoforms for a given gene on a given chromosome
    write_csv(save_as)
}

# ---- Codon counts ----
analysis_cds_codon_counts <- function(given   = 'data/CDS/Hsapiens-UCSC-hg38-validated.csv',
                                      save_as = 'data/CDS/Hsapiens-UCSC-hg38-validated-codon-counts.csv') {
  
  genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  
  read_csv(given, col_types = cols()) %>%
    # Extract and construct coding sequenes from CDS exon coordinates
    mutate(exon_seq = get_genomic_sequence(chr, strand, start, end, genome)) %>%
    arrange(gene, tx, exon) %>%
    group_by(gene, tx) %>%
    summarise(
      sequence   = str_c(exon_seq, collapse = ''),
      cds_length = str_length(sequence)
    ) %>%
    ungroup() %>%
    # Count the codons
    mutate(codon_counts = .count_codons(sequence)) %>%
    select(-sequence) %>% # remove for to reduce size
    unnest() %>% # expands the list of data frames in codon_counts to their own columns
    write_csv(save_as)
}


analysis_cds_nonsense_counts <-function(given   = "data/CDS/Hsapiens-UCSC-hg38-validated-codon-counts.csv",
                                        save_as = 'data/CDS/Hsapiens-UCSC-hg38-validated-nonsense-counts.csv'){ 
  
  counts %>%
    select(tx, gene, cds_length, codon, codon_count) %>%
    spread(codon, codon_count, fill = 0) %>%
    mutate(
      AG = 0,
      GA = CAG + CAA + CGA + (2*TGG), 
      GT = TCG + GAG + TAC + GAA + TCA + TGC + GGA,
      AC = TAT + TTA,
      GC = TAC + TCA,
      AT = AAG + TTG + TAT + AAA + TTA + AGA + TGT,
      P.AG = AG / (n_a + n_t),
      P.GA = GA / (n_g + n_c),
      P.GT = GT / (n_g + n_c),
      P.AC = AC / (n_a + n_t),
      P.GC = GC / (n_g + n_c),
      P.AT = AT / (n_a + n_t)
      
      # call from cosmiic the number of CT, GC ETC  to get the number of "trials" * the prob it's a nonesense 
      # to get the expected number of that kind of mutations causing nonsense. sum accross all mutations types
      # to get the total number of nonsense expected for that gene.
    ) %>%
    select(tx, gene, cds_length, AG, GA, GT, AC, GC, AT)
}





# ---- CDS analysis helper functions ----

# Given vector of coding sequences return a list of codon count data frames
.count_codons <- function(x, type = Biostrings::DNAString) {
  map(x, function(x) {
    Biostrings::codons(type(x)) %>%
      as.data.frame() %>%
      count(x) %>%
      rename(codon = x, codon_count = n) %>%
      mutate(codon_proportion = codon_count / sum(codon_count))
  })
}

# Given a CSV of CDS exon coordinates, check for expected features and save as CSV
.validate_cds <- function(cds, genome) {
  
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
      near(1, .n_stop_codons(dna)),      # only one stop codon
      str_detect(dna, 'TGA$|TAA$|TAG$')  # stop codon is at end
    )
  
  filter(cds, tx %in% valid$tx) 
}

.n_stop_codons <- function(x) {
  str_locate_all(x, 'TGA|TAA|TAG') %>% map_int(~sum(.[, 'start'] %% 3 == 1L))
}
