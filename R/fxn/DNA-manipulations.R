fxn$reading_frame <- function(x, offset = 1L) ((x - 1L) %% 3L) + 1L

fxn$dna_reverse_complement <- function(x, complement = c('ATGCYRSWKMBDHVN' = 'TACGRYSWMKVHDBN')) {
  from <- paste0(names(complement), str_to_lower(names(complement)))
  to   <- paste0(complement, str_to_lower(complement))
  chartr(from, to, stringi::stri_reverse(x))
}

fxn$dna_translate <- function(x) {
  CODE <- structure(c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", 
              "*", "*", "C", "C", "*", "W", "L", "L", "L", "L", "P", "P", "P", 
              "P", "H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M", 
              "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "V", 
              "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", 
              "G", "G"), .Names = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", 
                                    "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", 
                                    "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", 
                                    "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", 
                                    "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", 
                                    "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", 
                                    "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", 
                                    "GGT", "GGC", "GGA", "GGG"), alt_init_codons = character(0))
  
  Biostrings::DNAStringSet(x) %>%
    Biostrings::translate(if.fuzzy.codon = 'solve', genetic.code = CODE) %>%
    as.character()
}

fxn$extract_guide <- function(coding_context, target_position = 31, base_edit = 'c', spacing = 12, pam = '.G[GA]') {
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
