figures_messense_silent_nonsense <- function(given_veps = c(clinvar = 'data/ClinVar/VEP.csv',
                                                            cosmic  = 'data/COSMIC/VEP.csv',
                                                            exac    = 'data/ExAc/VEP.csv'),
                                             save_as = 'figures/messense-silent-nonsense.pdf') {
  
  type_order <- c('Nonsense', 'Missense', 'Silent')
  
  columns <- cols_only(
    chr           = col_character(), 
    gene          = col_character(), 
    exon_strand   = col_character(),
    variant_type  = col_character(), 
    splicing_type = col_character(),
    mutation_type = col_character(),
    ref_aa        = col_character(), 
    alt_aa        = col_character(),
    ref_cds       = col_character(), 
    alt_cds       = col_character(), 
    POS           = col_integer(), 
    REF           = col_character(), 
    ALT           = col_character()
  )
  
  clin <- read_csv(given_veps['clinvar'], col_types = columns) %>% distinct()
  cosm <- read_csv(given_veps['cosmic'],  col_types = columns) %>% distinct()
  exac <- read_csv(given_veps['exac'],    col_types = columns) %>% distinct()
  
  data <- bind_rows(lst(clin, cosm, exac), .id = 'dataset')
  
  summary <-
    data %>%
    filter(mutation_type %in% type_order) %>%
    count(dataset, mutation_type) %>%
    group_by(dataset) %>%
    mutate(
      total = sum(n),
      percent = n / total
    ) %>%
    ungroup()
  
  summary %>%
    mutate(
      mutation_type = ordered(mutation_type, levels = type_order),
      dataset = ordered(dataset, levels = c('exac', 'cosm', 'clin'), labels = c('ExAc', 'COSMIC', 'ClinVar'))
    ) %>%
    ggplot(aes(x = dataset, y = percent, fill = mutation_type)) +
    geom_col(alpha = 0.5, color = 'black') +
    geom_text(aes(label = scales::percent(round(percent, 2))), position = position_stack(vjust = 0.5)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c('red', 'lightblue', 'grey')) +
    coord_cartesian(ylim = c(-0.01, 1.01), expand = FALSE) +
    labs(x = '', y = '% of coding variants', fill = '') +
    theme_bw() +
    theme(aspect.ratio = 4, panel.border = element_blank(), panel.grid.major.x = element_blank()) +
    ggsave(save_as, width = 8, height = 10)
}
