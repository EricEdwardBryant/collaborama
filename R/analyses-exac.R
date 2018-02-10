analyses_exac_snv_counts <- function(given = 'data/ExAc/ExAC.r1.sites.vep.vcf.gz',
                                     save_as = 'data/ExAc/exac-snv-counts.csv') {
  
  callback <- function(x, pos) {
    x %>%
      # Only consider "passing" variants
      filter(FILTER == 'PASS') %>%
      mutate(
        mutation_class = case_when(
          (REF == 'A' & ALT == 'C') | (REF == 'T' & ALT == 'G') ~ 'A to C (T to G)',
          (REF == 'A' & ALT == 'G') | (REF == 'T' & ALT == 'C') ~ 'A to G (T to C)', # Transition
          (REF == 'A' & ALT == 'T') | (REF == 'T' & ALT == 'A') ~ 'A to T (T to A)',
          (REF == 'C' & ALT == 'A') | (REF == 'G' & ALT == 'T') ~ 'C to A (G to T)',
          (REF == 'C' & ALT == 'G') | (REF == 'G' & ALT == 'C') ~ 'C to G (G to C)',
          (REF == 'C' & ALT == 'T') | (REF == 'G' & ALT == 'A') ~ 'C to T (G to A)', # Transition
          TRUE ~ 'Other'
        )
      ) %>%
      count(mutation_class)
  }
  
  # Sum the count of each chunk
  exac <- 
    read_vcf(given, callback = callback) %>%
    group_by(mutation_class) %>%
    summarise(n = sum(n))

  write_csv(exac, save_as)
}

analyses_exac_vep <- function(given_vcf   = 'data/ExAc/ExAC.r1.sites.vep.vcf.gz',
                              given_cds   = 'data/CDS/Hsapiens-UCSC-hg19-validated.csv',
                              genome  = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                              save_as = 'data/ExAc/VEP.csv') {
  # Initialize output file
  header <- 'ID,chr,gene,exon_strand,variant_type,splicing_type,mutation_type,ref_aa,alt_aa,ref_cds,alt_cds,POS,REF,ALT,inframe_ref_start,inframe_ref_end,inframe_alt_start,inframe_alt_end'
  write_lines(header, save_as)
  
  cds <- read_csv(given_cds)
  
  vep_callback <- function(x, pos) {
    vcf <- x %>%
      filter(FILTER == 'PASS') %>% # Only worry about passing variants for now
      mutate(CHROM = str_c('chr', str_replace(CHROM, 'MT', 'M')))
    
    try(
      mutagenesis::predict_variant_effect(cds, vcf, genome) %>%
      select(
        ID, chr, gene, exon_strand, variant_type, splicing_type,
        mutation_type, ref_aa, alt_aa, ref_cds, alt_cds, POS, REF, ALT,
        inframe_ref_start, inframe_ref_end, inframe_alt_start, inframe_alt_end
      ) %>%
      distinct() %>%
      write_csv(save_as, append = TRUE) # keep adding to the file
    )
    
    # Every successful block
    return(data_frame(REF = 'A', ALT = 'G'))
  }
  
  read_vcf(given_vcf, callback = vep_callback)
  return(invisible())
}
