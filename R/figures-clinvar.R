figures_clinvar_summary <- function(given = c(all = 'data/ClinVar/ClinVar-summary-variant-types.csv',
                                              SNV = 'data/ClinVar/ClinVar-creating-vs-reverting.csv'),
                                    save_as = 'figures/Clinvar-summary.pdf') {

  # Color and ordering
  order = c(
    'SNV' = '#eab65e',
    'Insertion/Deletion' = '#d8d9de',
    'Duplication' = '#707887',
    # Transitions
    'C to T (G to A)' = '#b3a0c4', 
    'A to G (T to C)' = '#604488',
    # Transversions
    'C to G (G to C)' = '#d2e9c8',
    'C to A (G to T)' = '#b3dbb8',
    'A to C (T to G)' = '#87c17e',
    'A to T (T to A)' = '#60a865'
  )
  
  # Dataset
  all <- 
    read_csv(given['all'], col_types = cols()) %>%
    mutate(set = 'All variants') %>%
    filter(variant != 'Other')
  snv <- 
    read_csv(given['SNV'], col_types = cols()) %>% 
    rename(variant = group) %>%
    filter(variant != 'Misannotation', variant != 'Missing or ambiguous allele information') %>%
    mutate(set = 'SNV only')
  data <-
    bind_rows(all, snv) %>%
    mutate(variant = ordered(variant, levels = names(order)))

  # Figure
  ggplot(data, aes(y = percent / 100, x = set, fill = variant, label = variant)) +
    facet_grid(~set, scales = 'free_x') +
    geom_col(color = 'black', position = position_stack(), alpha = 0.5) + 
    geom_text(color = 'black', size = 5, fontface = 'bold', position = position_stack(vjust = 0.5)) +
    scale_color_manual(values = order) + 
    scale_fill_manual(values = order) + 
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    theme_bw() +
    guides(color = 'none', fill = 'none', label = 'none') +
    theme(
      panel.border = element_blank(), 
      panel.grid.major.x = element_blank(), 
      axis.text.y = element_text(size = 12), 
      strip.text = element_text(size = 12, face = 'bold'), 
      axis.title = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.text.x = element_blank(), 
      strip.background = element_blank(), 
      panel.spacing.x = unit(1, 'cm')
    ) +
    ggsave(save_as, width = 5, height = 8)
}


figures_clinvar_snv_counts <- function(given = 'data/ClinVar/ClinVar-creating-vs-reverting.csv') {
  read_csv(given)
}
