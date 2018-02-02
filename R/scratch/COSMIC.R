data_import$Cosmic_clean <- function(given_cosmic = 'data/COSMIC/CosmicMutantExport.tsv.gz', 
                                     given_cds = 'data/CDS/Hsapiens-UCSC-hg38.csv',
                                     save_as = 'data/COSMIC/Cosmic.csv.gz') {
  
  CDS_contiguous_intervals <- read_csv(given_cds, col_types = cols())
  
  given_cosmic %>%
    read_tsv(
      col_types = cols_only(
        #'ID_sample' = col_integer(),
        #'ID_tumour'  = col_integer(),
        'Primary site' = col_character(),
        'Mutation ID' = col_character(),
        'Primary histology' = col_character(),
        'Mutation zygosity' = col_character(),
        'Mutation Description' = col_character(),
        'Mutation CDS' = col_character(),
        'Mutation AA' = col_character(),
        'GRCh' = col_integer(),
        'Mutation genome position' = col_character(),
        'Mutation strand' = col_character()
        #'FATHMM prediction' = col_character(),
        #'FATHMM score' = col_double(),
        #'Mutation somatic status' = col_character()
      )
    ) %>%
    # Use hg38 coordinates and only consider frameshift (Indel) or substitutions
    filter(GRCh == 38) %>%
    # Extract chromosome coordinates from 'Mutation genome position' string
    mutate(
      chr    = str_c('chr', str_extract(`Mutation genome position`, '.*(?=:)')),
      strand = `Mutation strand`,
      start  = str_extract(`Mutation genome position`, '(?<=:).*(?=-)') %>% as.integer(),
      end    = str_extract(`Mutation genome position`, '(?<=-).*') %>% as.integer()
    ) %>%
    # Add gene name based on UCSC coordinates of CDS
    fuzzyjoin::genome_left_join(CDS_contiguous_intervals, by = c('chr', 'start', 'end')) %>%
    filter(strand.x == strand.y) %>% # Strand information should correspond to the gene's strand
    distinct() %>%
    # Define approximate cancer types and mutation classes
    mutate(
      cancer_type = case_when(
        str_detect(`Primary histology`, 'melanoma')                      ~ 'Malignant melanoma',
        `Primary site` == 'large_intestine'                              ~ 'Colorectal',
        `Primary site` == 'endometrium'                                  ~ 'Endometrial',
        `Primary site` == 'lung'                                         ~ 'Lung',
        `Primary site` == 'liver'                                        ~ 'Liver',
        `Primary site` == 'skin'                                         ~ 'Non-melanoma skin',
        `Primary site` == 'breast'                                       ~ 'Breast',
        `Primary site` == 'stomach'                                      ~ 'Stomach',
        `Primary site` %in% c('upper_aerodigestive_tract', 'oesophagus') ~ 'Upper aerodigestive',
        `Primary site` == 'haematopoietic_and_lymphoid_tissue'           ~ 'Blood',
        `Primary site` == 'prostate'                                     ~ 'Prostate',
        `Primary site` == 'pancreas'                                     ~ 'Pancreatic',
        `Primary site` == 'urinary_tract'                                ~ 'Bladder',
        `Primary site` == 'kidney'                                       ~ 'Kidney',
        `Primary histology` == 'glioma'                                  ~ 'Glioma',
        `Primary site` == 'ovary'                                        ~ 'Ovarian',
        `Primary site` == 'cervix'                                       ~ 'Cervical',
        `Primary site` == 'thyroid'                                      ~ 'Thyroid',
        `Primary site` == 'bone'                                         ~ 'Bone',
        TRUE ~ 'Other'
      ),
      mutation_class = case_when(
        str_detect(`Mutation Description`, 'Frameshift') ~ 'Frameshift',
        str_detect(`Mutation Description`, 'Missense')   ~ 'Missense',
        str_detect(`Mutation Description`, 'Nonsense')   ~ 'Nonsense',
        str_detect(`Mutation Description`, 'silent')     ~ 'Silent',
        `Mutation Description` == 'Nonstop extension'    ~ 'Nonstop extension',
        TRUE ~ 'Other'
      )
    ) %>%
    # Re-name and re-order columns
    select(
      gene, cancer_type, mutation_id = `Mutation ID`, mutation_class,
      chr = chr.x, strand = strand.x, start = start.x, end = end.x, 
      description = `Mutation Description`, mutation_cds = `Mutation CDS`, 
      mutation_aa = `Mutation AA`
    ) %>%
    distinct() %>%
    mutate(
      nt_from = mutation_cds %>% str_replace('^c.*[0-9]', '') %>% str_replace('>.*', ''),
      nt_to   = str_extract(mutation_cds, '(?<=([>])).*'),
      aa_from = str_extract(mutation_aa,  '(?<=(p[.])).*?(?=([0-9]))'),
      aa_to   = str_replace(mutation_aa,  '^p.*[0-9]', '')
    ) %>%
    write_csv(save_as)
}


data_import$Cosmic_mutation_subset <- function(given   = 'data/COSMIC/Cosmic.csv.gz',
                                               class   = 'Nonsense',
                                               save_as = 'data/COSMIC/Cosmic-nonsense.csv') {
  read_csv(given, col_types = cols()) %>%
    filter(mutation_class %in% class) %>%
    select(-mutation_class) %>%
    write_csv(save_as)
}