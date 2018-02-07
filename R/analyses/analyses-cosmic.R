analyses$cosmic_vep <- function(given_vcf  = 'data/COSMIC/CosmicCodingMuts.vcf.gz',
                                given_cds = 'data/CDS/Hsapiens-UCSC-hg38-validated.csv',
                                save_as   = 'data/COSMIC/VEP.csv',
                                genome    = BSgenome.Hsapiens.UCSC.hg38::Hsapiens) {
  cds <- read_csv(given_cds)
  
  # Harmonize chromosome names in VCF
  vcf <- read_vcf(given_vcf) %>%
    mutate(CHROM = str_c('chr', str_replace(CHROM, 'MT', 'M')))
  
  mutagenesis::predict_variant_effect(cds, vcf, genome) %>%
    write_csv(save_as)
}

analyses$COSMIC_vcf_to_csv <- function(given   = 'data/COSMIC/CosmicCodingMuts.vcf.gz',
                                       save_as = 'data/COSMIC/CosmicCodingMuts.csv') {
  fxn$read_vcf(given) %>%
    mutate(
      # Fix chromosome names to match names(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
      CHROM = str_c('chr', CHROM) %>% str_replace('MT', 'M'), 
      # Extract SNP flag from INFO column
      SNP   = str_detect(INFO, 'SNP;'), 
      INFO  = str_replace(INFO, 'SNP;', '')
    ) %>%
    # Parse INFO column into gene, strand, CDS, AA and count columns
    separate(INFO, into = c('gene', 'strand', 'CDS', 'AA', 'count'), sep = ';') %>%
    mutate_at(
      # Remove 'key=' and do type conversion if possible
      c('gene', 'strand', 'CDS', 'AA', 'count'), 
      funs(type.convert(str_replace(., pattern = '^.*=', replacement = ''), as.is = T))
    ) %>%
    write_csv(save_as)
}

analyses$COSMIC_nonsense <- function(given   = 'data/COSMIC/CosmicCodingMuts.csv',
                                     save_as = 'data/COSMIC/CosmicNonsense.csv') {
  # Nonsnse mutations will be annotated with a * at the end of the AA field 
  read_csv(given, col_types = cols()) %>%
    filter(str_detect(AA, '\\*$')) %>%
    write_csv(save_as)
}
