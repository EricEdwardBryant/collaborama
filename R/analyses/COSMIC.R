#analyses$COSMIC_DDR_candidates <- function(given = 'data/COSMIC/Nonsense.csv',
#                                           save_as) {
#  read_csv(given, col_types = cols()) %>%
#    filter(
#      C12PAM_nC == 1L | C13PAM_nC == 1L | C14PAM_nC == 1L | C15PAM_nC == 1L | C16PAM_nC == 1L
#    ) %>%
#    filter(!is.na(CtoT), !is.na(AtoG)) %>%
#
#    filter(guides_1C_noG, gene %in% DDR$gene) %>% 
#      mutate(searched = ref_context) %>% 
#      add_RFLP(width = 30) %>% 
#      add_RFLP(width = 30, recognizes = 't') %>%
#      select(-searched)
#    
#    write_csv(results, save_as)
#}

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
