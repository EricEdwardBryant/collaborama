data_import$snv_to_stop <- function() {
  code <- 
    Biostrings::GENETIC_CODE %>%
    enframe(name = 'codon', value = 'aa') #%>%
    #mutate(codon_rev = fxn$dna_reverse_complement(codon_fwd))
  
  mutation_types <- tribble(
    ~from, ~to, ~mutation_symbol, ~mutation_type,
    'A', 'T'
  )
  
  snv <- code %>%
    mutate(
      snv_T1 = str_replace(codon, '^.',       'T'),
      snv_T2 = str_replace(codon, '(?<=.).',  'T'),
      snv_T3 = str_replace(codon, '.$',       'T'),
      snv_A1 = str_replace(codon, '^.',       'A'),
      snv_A2 = str_replace(codon, '(?<=.).',  'A'),
      snv_A3 = str_replace(codon, '.$',       'A'),
      snv_G1 = str_replace(codon, '^.',       'G'),
      snv_G2 = str_replace(codon, '(?<=.).',  'G'),
      snv_G3 = str_replace(codon, '.$',       'G'),
      snv_C1 = str_replace(codon, '^.',       'C'),
      snv_C2 = str_replace(codon, '(?<=.).',  'C'),
      snv_C3 = str_replace(codon, '.$',       'C')
    ) %>%
    gather('snv_type', 'snv_codon', starts_with('snv')) %>%
    mutate(
      from = case_when(
        str_detect(snv_type, '1') ~ str_sub(codon, 1, 1),
        str_detect(snv_type, '2') ~ str_sub(codon, 2, 2),
        str_detect(snv_type, '3') ~ str_sub(codon, 3, 3)
      ),
      to = str_sub(snv_type, 5, 5),
      snv_aa = fxn$dna_translate(snv_codon)
    ) %>%
    filter(aa != snv_aa)
  
  
  
  to_stop   <- filter(snv, snv_aa == '*')
  from_stop <- filter(snv, aa == '*')

}
