data_import$ClinVar_vcf_to_csv <- function(given   = 'data/ClinVar/clinvar.vcf.gz',
                                           save_as = 'data/COSMIC/ClinVar.csv') {
  
  read_vcf(given) %>%
    mutate(
      # Fix chromosome names to match names(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
      CHROM        = str_c('chr', CHROM) %>% str_replace('MT', 'M'),
      # Extract Info
      AF_ESP       = str_extract(INFO, '(?<=(AF_ESP=)).*?(?=;)'),
      AF_EXAC      = str_extract(INFO, '(?<=(AF_EXAC=)).*?(?=;)'),
      AF_TGP       = str_extract(INFO, '(?<=(AF_TGP=)).*?(?=;)'),
      ALLELEID     = str_extract(INFO, '(?<=(ALLELEID=)).*?(?=;)'),
      CLNDN        = str_extract(INFO, '(?<=(CLNDN=)).*?(?=;)'),
      CLNDNINCL    = str_extract(INFO, '(?<=(CLNDNINCL=)).*?(?=;)'),
      CLNDISDB     = str_extract(INFO, '(?<=(CLNDISDB=)).*?(?=;)'),
      CLNDISDBINCL = str_extract(INFO, '(?<=(CLNDISDBINCL=)).*?(?=;)'),
      CLNHGVS      = str_extract(INFO, '(?<=(CLNHGVS=)).*?(?=;)'),
      CLNREVSTAT   = str_extract(INFO, '(?<=(CLNREVSTAT=)).*?(?=;)'),
      CLNSIG       = str_extract(INFO, '(?<=(CLNSIG=)).*?(?=;)'),
      CLNSIGINCL   = str_extract(INFO, '(?<=(CLNSIGINCL=)).*?(?=;)'),
      CLNVC        = str_extract(INFO, '(?<=(CLNVC=)).*?(?=;)'),
      CLNVCSO      = str_extract(INFO, '(?<=(CLNVCSO=)).*?(?=;)'),
      CLNVI        = str_extract(INFO, '(?<=(CLNVI=)).*?(?=;)'),
      DBVARID      = str_extract(INFO, '(?<=(DBVARID=)).*?(?=;)'),
      GENEINFO     = str_extract(INFO, '(?<=(GENEINFO=)).*?(?=;)'),
      MC           = str_extract(INFO, '(?<=(MC=)).*?(?=;)'),
      ORIGIN       = str_extract(INFO, '(?<=(ORIGIN=)).*?(?=;)'),
      RS           = str_extract(INFO, '(?<=(RS=)).*?(?=;)'),
      SSR          = str_extract(INFO, '(?<=(SSR=)).*?(?=;)')
    ) %>%
    select(-INFO) %>%
    mutate_if(is.character, type.convert, as.is = TRUE) %>%
    write_csv(save_as)
}

data_import$ClinVar_nonsense <- function(given   = 'data/ClinVar/ClinVar.csv',
                                         save_as = 'data/ClinVar/ClinVar-Nonsense.csv') {
  # Nonsnse mutations will be annotated with a * at the end of the AA field 
  read_csv(given, col_types = cols()) %>%
    filter(str_detect(MC, 'nonsense')) %>%
    write_csv(save_as)
}
