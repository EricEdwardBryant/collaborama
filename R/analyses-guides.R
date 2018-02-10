analysis_design_guides_to_model_transitions_with_BE3_and_ABE <- function(given   = 'data/COSMIC/Cosmic-Coding.csv',
                                                                         save_as = 'data/Guides/Cosmic-Coding-BE3-ABE.csv') {
  read_csv(given, col_types = cols()) %>%
    .design_guides_to_model_transitions_with_BE3_and_ABE() %>%
    write_csv(save_as)
}

# Design guides to model C>T variants using BE3 and ABE
#
# @param given Path to a CSV file with columns [CHROM, POS, REF, ALT]
# 
# @Details Only variants that cause SNV transitions are currently considered

.design_guides_to_model_transitions_with_BE3_and_ABE <- function(df,
                                                                 # Order of spacing used to prioritize final guides
                                                                 ABE_spacing = c(14, 13, 15, 12, 16),
                                                                 BE3_spacing = c(14, 13, 15, 12, 16),
                                                                 BE3_PAM_pattern = c('.GG', '.GA', '.GCG', '...[AG][AG]T'), # full PAM diversity captured
                                                                 ABE_PAM_pattern = '.GG',
                                                                 flanking = 30,
                                                                 Genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens) {
  message('Extracting genomic sequence context...')
  result <- 
    as_tibble(df) %>%
    filter(
      (REF == 'C' & ALT == 'T') | 
        (REF == 'G' & ALT == 'A') |
        (REF == 'T' & ALT == 'C') |
        (REF == 'A' & ALT == 'G')
    ) %>% 
    mutate(
      # Last coordinate of REF (only relevant if REF is wider than one base)
      END = POS + str_length(REF) - 1L,
      
      # Get genomic sequence for 30 bp upstream and downstream of REF. Reference and alternative sequence in lowercase
      five_prime  = iSTOP:::get_genomic_sequence(POS, add_5prime = flanking, add_3prime = -1,       Genome, CHROM, strand = '+'),
      three_prime = iSTOP:::get_genomic_sequence(END, add_5prime = -1,       add_3prime = flanking, Genome, CHROM, strand = '+'),
      ref_check   = iSTOP:::get_genomic_sequence(POS, add_5prime =  0,       add_3prime = 0,        Genome, CHROM, strand = '+'),
      ref_plus    = str_c(five_prime, str_to_lower(REF), three_prime),
      alt_plus    = str_c(five_prime, str_to_lower(ALT), three_prime),
      ref_minus   = reverse_complement(ref_plus),
      alt_minus   = reverse_complement(alt_plus),
      
      # Annotate Gc status
      Gc = str_detect(ref_plus, 'Gc') | str_detect(ref_minus, 'Gc'),
      
      # v--- NOT NECESSARY ---v    
      # Determine the frame of the first reference position
      ref_frame    = ((POS - 1L) %% 3L) + 1L,
      # Translate all three frames of the alternative sequence context
      frame1_plus  = translate(str_sub(alt_plus,  1L, 57L)),
      frame2_plus  = translate(str_sub(alt_plus,  2L, 58L)),
      frame3_plus  = translate(str_sub(alt_plus,  3L, 59L)),
      frame1_minus = translate(str_sub(alt_minus, 1L, 57L)),
      frame2_minus = translate(str_sub(alt_minus, 2L, 58L)),
      frame3_minus = translate(str_sub(alt_minus, 3L, 59L)),
      # Determine which frame had the expected nonsense substitution 
      # (just an estimate, will break if splice junction occurs in codon)
      frame_plus = case_when(
        str_detect(frame1_plus, '^.{9,10}\\*') ~ 1L,
        str_detect(frame2_plus, '^.{9,10}\\*') ~ 2L,
        str_detect(frame3_plus, '^.{9,10}\\*') ~ 3L
      ),
      frame_minus = case_when(
        str_detect(frame1_minus, '^.{9,10}\\*') ~ 1L,
        str_detect(frame2_minus, '^.{9,10}\\*') ~ 2L,
        str_detect(frame3_minus, '^.{9,10}\\*') ~ 3L
      )
      # ^--- NOT NECESSARY ---^
    )
  
  # Constuct all BE3 editing guides using reference sequence context (NGG/NGA)
  # Count the number of Cs in BE3 editing window without an upstream G
  message('Designing BE3 guides...')
  BE3_model_cols  <- str_c('BE3_model_',  BE3_spacing)  # order used for coalescing guides
  BE3_revert_cols <- str_c('BE3_revert_',  BE3_spacing)  # order used for coalescing guides
  for (spacing in BE3_spacing) {
    BE3_model  <- str_c('BE3_model_', spacing)
    BE3_revert <- str_c('BE3_revert_', spacing)
    nC_model   <- str_c(BE3_model, '_nC')
    nC_revert  <- str_c(BE3_revert, '_nC')
    result[[BE3_model]]  <- fxn$extract_guide(result$ref_plus, target_position = flanking + 1, base_edit = 'c', spacing = spacing, pam = BE3_PAM_pattern)
    result[[BE3_revert]] <- fxn$extract_guide(result$alt_plus, target_position = flanking + 1, base_edit = 'c', spacing = spacing, pam = BE3_PAM_pattern)
    # Count the number of Cs in BE3 editing window with no upstream G (inhibitory to editing)
    result[[nC_model]]  <- str_sub(result[[BE3_model]],  3, 8) %>% str_to_upper() %>% str_count('(?<=([^G]))C') # C preceded by !G
    result[[nC_revert]] <- str_sub(result[[BE3_revert]], 3, 8) %>% str_to_upper() %>% str_count('(?<=([^G]))C') # C preceded by !G
    
  }
  
  # Constuct all ABE editing guides using alternate sequence context (NGG)
  # Count the number of As in ABE editing window
  message('Designing ABE guides...')
  ABE_model_cols  <- str_c('ABE_model_', ABE_spacing)  # order used for coalescing guides
  ABE_revert_cols <- str_c('ABE_revert_', ABE_spacing)  # order used for coalescing guides
  for (spacing in ABE_spacing) {
    ABE_model  <- str_c('ABE_model_', spacing)
    ABE_revert <- str_c('ABE_revert_', spacing)
    nA_model   <- str_c(ABE_model, '_nA')
    nA_revert  <- str_c(ABE_revert, '_nA')
    result[[ABE_model]]  <- .extract_guide(result$ref_plus, target_position = flanking + 1, base_edit = 'a', spacing = spacing, pam = ABE_PAM_pattern)
    result[[ABE_revert]] <- .extract_guide(result$alt_plus, target_position = flanking + 1, base_edit = 'a', spacing = spacing, pam = ABE_PAM_pattern)
    # Count the number of As in the ABE editing window
    result[[nA_model]]   <- str_sub(result[[ABE_model]],  4, 8) %>% str_to_upper() %>% str_count('A')
    result[[nA_revert]]  <- str_sub(result[[ABE_revert]], 4, 8) %>% str_to_upper() %>% str_count('A')
  }
  
  # Coalesce counts and guides
  message('Selecting preferred guides for BE3 and ABE...')
  result$BE3_model  <- coalesce(!!! result[, BE3_model_cols])
  result$BE3_revert <- coalesce(!!! result[, BE3_revert_cols])
  result$ABE_model  <- coalesce(!!! result[, ABE_model_cols])
  result$ABE_revert <- coalesce(!!! result[, ABE_revert_cols])
  
  
  # Safe BE3 guides are those that have only one C (excluding Gc)
  # The Gc gets added as it was not counted previously
  result$BE3_model_safe <- coalesce(!!! lapply(BE3_model_cols, function(x) {
    ifelse(result[[str_c(x, '_nC')]] + result$Gc == 1L, result[[x]], NA)
  }))
  result$BE3_revert_safe <- coalesce(!!! lapply(BE3_revert_cols, function(x) {
    ifelse(result[[str_c(x, '_nC')]] + result$Gc == 1L, result[[x]], NA)
  }))
  
  # Safe ABE guides are those that have only one A
  result$ABE_model_safe <- coalesce(!!! lapply(ABE_model_cols, function(x) {
    ifelse(result[[str_c(x, '_nA')]] == 1L, result[[x]], NA)
  }))
  result$ABE_revert_safe <- coalesce(!!! lapply(ABE_revert_cols, function(x) {
    ifelse(result[[str_c(x, '_nA')]] == 1L, result[[x]], NA)
  }))
  
  # Return result without flanking sequences used to construct ref and alt context
  select(result, -five_prime, -three_prime)
}


.extract_guide <- function(coding_context, target_position = 31, base_edit = 'c', spacing = 12, pam = '.G[GA]') {
  guide_lead <- 19 - spacing
  contx_lead <- 30 - guide_lead
  result <- rep(NA_character_, length(coding_context))
  for (ipam in pam) {
    pattern    <- str_c('(?<=(.{', contx_lead, '})).{', guide_lead, '}', base_edit, '.{', spacing, '}', ipam)
    coding     <- str_extract(coding_context, pattern)
    non_coding <- str_extract(fxn$dna_reverse_complement(coding_context), pattern)
    result     <- coalesce(result, coding, non_coding)
  }
  return(result)
}
