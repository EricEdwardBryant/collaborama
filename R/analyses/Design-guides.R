analyses$Design_guides_to_model_transitions_with_BE3_and_ABE <- function(given   = 'data/COSMIC/Cosmic-Coding.csv',
                                                                         save_as = 'data/Guides/Cosmic-Coding-BE3-ABE.csv') {
  read_csv(given, col_types = cols()) %>%
    fxn$Design_guides_to_model_transitions_with_BE3_and_ABE() %>%
    write_csv(save_as)
}
