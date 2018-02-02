# Question:
# What is the general distribution of variant types in ClinVar?
analyses$clinvar_summary_of_variant_types <- function(given = 'data/ClinVar/ClinVar.csv',
                                                      save_as = 'data/ClinVar/ClinVar-summary-variant-types.csv') {
  read_csv(given, col_types = cols()) %>%
    mutate(
      variant = case_when(
        CLNVC == 'single_nucleotide_variant'                     ~ 'SNV',
        CLNVC %in% c('Deletion', 'Insertion', 'Indel')           ~ 'Insertion/Deletion',
        CLNVC %in% c('Variation', 'Microsatellite', 'Inversion') ~ 'Other',
        TRUE ~ CLNVC
      )
    ) %>%
    group_by(variant) %>% 
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(percent = round((n / sum(n)) * 100, digits = 1)) %>%
    arrange(desc(percent)) %>%
    write_csv(save_as)
}

analyses$clinvar_vep <- function(given_vcf = 'data/ClinVar/clinvar.vcf.gz',
                                 given_cds = 'data/CDS/Hsapiens-UCSC-hg38-validated.csv',
                                 save_as   = 'data/ClinVar/VEP.csv',
                                 genome    = BSgenome.Hsapiens.UCSC.hg38::Hsapiens) {
  cds <- read_csv(given_cds)
  
  # Harmonize chromosome names in VCF
  vcf <- read_vcf(given_vcf) %>%
    mutate(CHROM = str_c('chr', str_replace(CHROM, 'MT', 'M')))

  mutagenesis::predict_variant_effect(cds, vcf, genome) %>%
    write_csv(save_as)
}

# Question:
# Among SNVs What is the general distribution of transitions versus transversions?
#
# Question:
# How many variants can be modeled or corrected by BE3 and how many by ABE?
#
# 1  A>C (T>G)
# 2  A>G (T>C)
# 3  A>T (T>A)
# 4  C>A (G>T)
# 5  C>G (G>C)
# 6  C>T (G>A)

analyses$clinvar_summary_of_creating_vs_reverting <- function(given   = 'data/ClinVar/ClinVar.csv',
                                                              save_as = 'data/ClinVar/ClinVar-creating-vs-reverting.csv',
                                                              Genome  = BSgenome.Hsapiens.UCSC.hg38::Hsapiens) {
    read_csv(given, col_types = cols()) %>%
      select(CHROM, POS, ID, REF, ALT, CLNVC) %>%
      filter(CLNVC == 'single_nucleotide_variant') %>% 
      mutate(
        group = case_when(
          is.na(ALT) | ALT == 'N'                               ~ 'Missing or ambiguous allele information',
          str_length(REF) != 1 | str_length(ALT) != 1           ~ 'Misannotation',
          (REF == 'A' & ALT == 'C') | (REF == 'T' & ALT == 'G') ~ 'A to C (T to G)',
          (REF == 'A' & ALT == 'G') | (REF == 'T' & ALT == 'C') ~ 'A to G (T to C)', # Transition
          (REF == 'A' & ALT == 'T') | (REF == 'T' & ALT == 'A') ~ 'A to T (T to A)',
          (REF == 'C' & ALT == 'A') | (REF == 'G' & ALT == 'T') ~ 'C to A (G to T)',
          (REF == 'C' & ALT == 'G') | (REF == 'G' & ALT == 'C') ~ 'C to G (G to C)',
          (REF == 'C' & ALT == 'T') | (REF == 'G' & ALT == 'A') ~ 'C to T (G to A)'  # Transition
        ),
        ts_tv = case_when(
          group %in% c('A to G (T to C)', 'C to T (G to A)')                                       ~ 'Transition',
          group %in% c('A to C (T to G)', 'A to T (T to A)', 'C to A (G to T)', 'C to G (G to C)') ~ 'Transversion'
        ),
        editor = case_when(group == 'A to G (T to C)' ~ 'ABE', group == 'C to T (G to A)' ~ 'BE3'),
      
        
        # Get genomic sequence for 30 bp upstream and downstream of REF. Reference and alternative sequence in lowercase
        five_prime  = get_genomic_sequence(CHROM, '+', start = POS - 30, end = POS - 1,  Genome),
        three_prime = get_genomic_sequence(CHROM, '+', start = POS + 1,  end = POS + 30, Genome),
        ref_check   = get_genomic_sequence(CHROM, '+', start = POS, end = POS, Genome),
        ref_plus    = str_c(five_prime, str_to_lower(REF), three_prime),
        alt_plus    = str_c(five_prime, str_to_lower(ALT), three_prime),
        ref_minus   = reverse_complement(ref_plus),
        alt_minus   = reverse_complement(alt_plus),
        
        # Determine availability of PAM for transition mutations
        Create_with_BE3_NGG    = ts_tv == 'Transition' & (str_detect(ref_plus, 'c.{12,16}.GG')          | str_detect(ref_minus, 'c.{12,16}.GG')),
        Create_with_BE3_NGA    = ts_tv == 'Transition' & (str_detect(ref_plus, 'c.{12,16}.GA')          | str_detect(ref_minus, 'c.{12,16}.GA')),
        Create_with_BE3_NGCG   = ts_tv == 'Transition' & (str_detect(ref_plus, 'c.{12,16}.GCG')         | str_detect(ref_minus, 'c.{12,16}.GCG')),
        Create_with_BE3_NGAG   = ts_tv == 'Transition' & (str_detect(ref_plus, 'c.{12,16}.GAG')         | str_detect(ref_minus, 'c.{12,16}.GAG')),
        Create_with_BE3_NNGRRT = ts_tv == 'Transition' & (str_detect(ref_plus, 'c.{12,16}..G[AG][AG]T') | str_detect(ref_minus, 'c.{12,16}..G[AG][AG]T')),
        Create_with_BE3_NNNRRT = ts_tv == 'Transition' & (str_detect(ref_plus, 'c.{12,16}...[AG][AG]T') | str_detect(ref_minus, 'c.{12,16}...[AG][AG]T')),
        Revert_with_BE3_NGG    = ts_tv == 'Transition' & (str_detect(alt_plus, 'c.{12,16}.GG')          | str_detect(alt_minus, 'c.{12,16}.GG')),
        Revert_with_BE3_NGA    = ts_tv == 'Transition' & (str_detect(alt_plus, 'c.{12,16}.GA')          | str_detect(alt_minus, 'c.{12,16}.GA')),
        Revert_with_BE3_NGCG   = ts_tv == 'Transition' & (str_detect(alt_plus, 'c.{12,16}.GCG')         | str_detect(alt_minus, 'c.{12,16}.GCG')),
        Revert_with_BE3_NGAG   = ts_tv == 'Transition' & (str_detect(alt_plus, 'c.{12,16}.GAG')         | str_detect(alt_minus, 'c.{12,16}.GAG')),
        Revert_with_BE3_NNGRRT = ts_tv == 'Transition' & (str_detect(alt_plus, 'c.{12,16}..G[AG][AG]T') | str_detect(alt_minus, 'c.{12,16}..G[AG][AG]T')),
        Revert_with_BE3_NNNRRT = ts_tv == 'Transition' & (str_detect(alt_plus, 'c.{12,16}...[AG][AG]T') | str_detect(alt_minus, 'c.{12,16}...[AG][AG]T')),
        Create_with_ABE_NGG    = ts_tv == 'Transition' & (str_detect(ref_plus, 'a.{12,16}.GG')          | str_detect(ref_minus, 'a.{12,16}.GG')),
        Create_with_ABE_NGA    = ts_tv == 'Transition' & (str_detect(ref_plus, 'a.{12,16}.GA')          | str_detect(ref_minus, 'a.{12,16}.GA')),
        Create_with_ABE_NGCG   = ts_tv == 'Transition' & (str_detect(ref_plus, 'a.{12,16}.GCG')         | str_detect(ref_minus, 'a.{12,16}.GCG')),
        Create_with_ABE_NGAG   = ts_tv == 'Transition' & (str_detect(ref_plus, 'a.{12,16}.GAG')         | str_detect(ref_minus, 'a.{12,16}.GAG')),
        Create_with_ABE_NNGRRT = ts_tv == 'Transition' & (str_detect(ref_plus, 'a.{12,16}..G[AG][AG]T') | str_detect(ref_minus, 'a.{12,16}..G[AG][AG]T')),
        Create_with_ABE_NNNRRT = ts_tv == 'Transition' & (str_detect(ref_plus, 'a.{12,16}...[AG][AG]T') | str_detect(ref_minus, 'a.{12,16}...[AG][AG]T')),
        Revert_with_ABE_NGG    = ts_tv == 'Transition' & (str_detect(alt_plus, 'a.{12,16}.GG')          | str_detect(alt_minus, 'a.{12,16}.GG')),
        Revert_with_ABE_NGA    = ts_tv == 'Transition' & (str_detect(alt_plus, 'a.{12,16}.GA')          | str_detect(alt_minus, 'a.{12,16}.GA')),
        Revert_with_ABE_NGCG   = ts_tv == 'Transition' & (str_detect(alt_plus, 'a.{12,16}.GCG')         | str_detect(alt_minus, 'a.{12,16}.GCG')),
        Revert_with_ABE_NGAG   = ts_tv == 'Transition' & (str_detect(alt_plus, 'a.{12,16}.GAG')         | str_detect(alt_minus, 'a.{12,16}.GAG')),
        Revert_with_ABE_NNGRRT = ts_tv == 'Transition' & (str_detect(alt_plus, 'a.{12,16}..G[AG][AG]T') | str_detect(alt_minus, 'a.{12,16}..G[AG][AG]T')),
        Revert_with_ABE_NNNRRT = ts_tv == 'Transition' & (str_detect(alt_plus, 'a.{12,16}...[AG][AG]T') | str_detect(alt_minus, 'a.{12,16}...[AG][AG]T')),
        Create_with_BE3 = Create_with_BE3_NGG | Create_with_BE3_NGA | Create_with_BE3_NGCG | Create_with_BE3_NGAG | Create_with_BE3_NNGRRT | Create_with_BE3_NNNRRT,
        Revert_with_BE3 = Revert_with_BE3_NGG | Revert_with_BE3_NGA | Revert_with_BE3_NGCG | Revert_with_BE3_NGAG | Revert_with_BE3_NNGRRT | Revert_with_BE3_NNNRRT,
        Create_with_ABE = Create_with_ABE_NGG | Create_with_ABE_NGA | Create_with_ABE_NGCG | Create_with_ABE_NGAG | Create_with_ABE_NNGRRT | Create_with_ABE_NNNRRT,
        Revert_with_ABE = Revert_with_ABE_NGG | Revert_with_ABE_NGA | Revert_with_ABE_NGCG | Revert_with_ABE_NGAG | Revert_with_ABE_NNGRRT | Revert_with_ABE_NNNRRT
      ) %>%
      group_by(group, ts_tv, editor) %>% 
      summarise(
        count                  = n(),
        Create_with_BE3        = sum(Create_with_BE3),
        Revert_with_BE3        = sum(Revert_with_BE3),
        Create_with_ABE        = sum(Create_with_ABE),
        Revert_with_ABE        = sum(Revert_with_ABE),
        Create_with_BE3_NGG    = sum(Create_with_BE3_NGG),
        Create_with_BE3_NGA    = sum(Create_with_BE3_NGA),
        Create_with_BE3_NGCG   = sum(Create_with_BE3_NGCG),
        Create_with_BE3_NGAG   = sum(Create_with_BE3_NGAG),
        Create_with_BE3_NNGRRT = sum(Create_with_BE3_NNGRRT),
        Create_with_BE3_NNNRRT = sum(Create_with_BE3_NNNRRT),
        Revert_with_BE3_NGG    = sum(Revert_with_BE3_NGG),
        Revert_with_BE3_NGA    = sum(Revert_with_BE3_NGA),
        Revert_with_BE3_NGCG   = sum(Revert_with_BE3_NGCG),
        Revert_with_BE3_NGAG   = sum(Revert_with_BE3_NGAG),
        Revert_with_BE3_NNGRRT = sum(Revert_with_BE3_NNGRRT),
        Revert_with_BE3_NNNRRT = sum(Revert_with_BE3_NNNRRT),
        Create_with_ABE_NGG    = sum(Create_with_ABE_NGG),
        Create_with_ABE_NGA    = sum(Create_with_ABE_NGA),
        Create_with_ABE_NGCG   = sum(Create_with_ABE_NGCG),
        Create_with_ABE_NGAG   = sum(Create_with_ABE_NGAG),
        Create_with_ABE_NNGRRT = sum(Create_with_ABE_NNGRRT),
        Create_with_ABE_NNNRRT = sum(Create_with_ABE_NNNRRT),
        Revert_with_ABE_NGG    = sum(Revert_with_ABE_NGG),
        Revert_with_ABE_NGA    = sum(Revert_with_ABE_NGA),
        Revert_with_ABE_NGCG   = sum(Revert_with_ABE_NGCG),
        Revert_with_ABE_NGAG   = sum(Revert_with_ABE_NGAG),
        Revert_with_ABE_NNGRRT = sum(Revert_with_ABE_NNGRRT),
        Revert_with_ABE_NNNRRT = sum(Revert_with_ABE_NNNRRT)
      ) %>%
      ungroup() %>%
      mutate(percent = round((count / sum(count)) * 100, digits = 1)) %>%
      select(group, ts_tv, editor, count, percent, everything()) %>%
      arrange(desc(percent)) %>%
      write_csv(save_as)
}
