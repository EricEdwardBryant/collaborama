data_import$CDS_hg38 <- function(save_as = 'data/CDS/Hsapiens-UCSC-hg38.csv') {
  if (!file.exists(save_as)) {
    iSTOP::CDS_Hsapiens_UCSC_hg38() %>% 
      write_csv(save_as)
  }
}

data_import$CDS_hg19 <- function(save_as = 'data/CDS/Hsapiens-UCSC-hg19.csv') {
  if (!file.exists(save_as)) {
    iSTOP::CDS(
      tx      = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz',
      gene    = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz',
      tx_cols = c(
        'tx', 'chr', 'strand', 'tx_start', 'tx_end', 'cds_start', 'cds_end',
        'exon_count', 'exon_start', 'exon_end', 'protein_id', 'align_id'),
      gene_cols = c(
        'tx', 'mRNA', 'spID', 'spDisplayID', 'gene', 'refseq', 'protAcc',
        'description', 'rfamAcc', 'tRnaName')
    ) %>% 
      write_csv(save_as)
  }
}

data_import$CDS_contiguous_intervals <- function(given = 'data/CDS/Hsapiens-UCSC-hg38.csv', 
                                                 save_as = 'data/CDS/Hsapiens-UCSC-hg38-contiguous-intervals.csv') {
  # All distinct CDS coordinates for genes
  CDS <- 
    read_csv(given, col_types = cols()) %>% 
    select(gene, chr, strand, start, end) %>% 
    distinct()
  
  # These are all the contiguous intervals of coding sequence for every gene
  IRanges::reduce(with(CDS, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), strand))) %>%
    # The rest is just cleaning up and re-annotating with gene names
    IRanges::as.data.frame() %>%
    mutate(chr = as.character(seqnames)) %>%
    select(chr, strand, start, end) %>%
    fuzzyjoin::genome_left_join(CDS, c('chr', 'start', 'end')) %>% 
    filter(strand.x == strand.y) %>%
    select(gene, chr = chr.x, strand = strand.x, start = start.x, end = end.x) %>%
    distinct() %>%
    write_csv(save_as)
}



