locate_PAM_generic <- function(coordinates,
                               genome,
                               spacing = iSTOP::PAM_spacing(),
                               PAM = iSTOP::PAM_patterns_default(),
                               flanking = 150,
                               keep_PAM = FALSE) {
  
  assertthat::assert_that(
    length(flanking) == 1,
    is.numeric(flanking),
    msg = 'Please ensure that `flanking` is a single number. For more help, see ?locate_PAM'
  )
  
  assertthat::assert_that(
    !is.null(names(PAM)),
    !any(names(PAM) == ''),
    all(map_lgl(PAM, ~all(utils::hasName(., c('pattern', 'width'))))),
    all(map_lgl(PAM, ~is.numeric(.$width))),
    all(map_lgl(PAM, ~is.character(.$pattern))),
    msg = 'Please ensure that all `PAM` patterns are named and have a single pattern and single numeric width. For more help, see ?locate_PAM'
  )
  
  assertthat::assert_that(
    all(utils::hasName(coordinates, c('genome_coord', 'chr', 'sg_strand'))),
    msg = 'Coordinates table must have columns named "genome_coord", "chr" and "sg_strand"'
  )
  
  if (nrow(coordinates) < 1) return(invisible(coordinates))
  
  targetable_coords    <- filter(coordinates, !is.na(genome_coord))
  
  sequences <-
    targetable_coords %>%
    mutate(
      # Get flanking genomic sequence context
      sequence = iSTOP:::get_genomic_sequence(
        at = genome_coord,
        add_5prime = flanking,
        add_3prime = flanking,
        genome,
        chr,
        sg_strand
      ),
      # make targeted base lower case
      sequence =  stringr::str_c(
        stringr::str_sub(sequence, end   = flanking),                                                     # LHS
        stringr::str_sub(sequence, start = flanking + 1, end = flanking + 1) %>% stringr::str_to_lower(), # C -> c
        stringr::str_sub(sequence, start = flanking + 2)                                                  # RHS
      )
    )
  
  # Add an sgSTOP column for each PAM
  for (i in 1:length(PAM)) {
    col_name_guide   <- names(PAM)[i]
    col_name_spacing <- names(PAM)[i] %>% stringr::str_c('_spacing')
    
    # Initialize columns
    sequences[[col_name_guide]]   <- NA_character_
    sequences[[col_name_spacing]] <- NA_character_
    for (j in 1:length(spacing)) {
      
      # Coalesce prioritizing previously detected guides
      guide_sequence <-
        coalesce(
          sequences[[col_name_guide]],
          iSTOP:::sgSTOP(sequences$sequence, pattern = PAM[[i]]$pattern, base_edit = '[atgc]', spacing = spacing[[j]], width = 20 + PAM[[i]]$width)
        )
      
      guide_spacing <-
        coalesce(
          sequences[[col_name_spacing]],
          ifelse(!is.na(guide_sequence), names(spacing)[j], NA_character_)
        )
      
      # Optionally remove PAM sequence
      if (!keep_PAM) guide_sequence <- stringr::str_sub(guide_sequence, start = 1L, end = 20L)
      
      # Update
      sequences[[col_name_guide]]   <- guide_sequence
      sequences[[col_name_spacing]] <- guide_spacing
    }
  }
  
  # If not all of the PAM columns are NA then there is at least 1 match
  sequences$match_any <- !apply(apply(sequences[, names(PAM)], 2, is.na), 1, all)
  
  return(sequences)
}
