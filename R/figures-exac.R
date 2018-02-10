figures_exac_snv_counts <- function(given = 'data/ExAc/exac-snv-counts.csv', 
                                    save_as = 'figures/exac-snv-counts.pdf') {
  # Color and ordering
  order = c(
    #'SNV' = '#eab65e',
    #'Insertion/Deletion' = '#d8d9de',
    #'Duplication' = '#707887',
    # Transitions
    'C to T (G to A)' = '#b3a0c4',
    'A to G (T to C)' = '#604488',
    # Transversions
    'C to G (G to C)' = '#d2e9c8',
    'C to A (G to T)' = '#b3dbb8',
    'A to C (T to G)' = '#87c17e',
    'A to T (T to A)' = '#60a865'
  )
  
  read_csv(given) %>% 
    filter(mutation_class != 'Other') %>% 
    mutate(
      set = 'SNV only',
      mutation_class = ordered(mutation_class, levels = names(order))
    ) %>%
    ggplot(aes(x = set, y = n / sum(n), fill = mutation_class, label = str_c(mutation_class, '\n', scales::percent(n / sum(n))))) + 
    facet_grid(~set, scales = 'free_x') +
    geom_col(color = 'black', position = position_stack(), alpha = 0.5) + 
    geom_text(color = 'black', size = 5, fontface = 'bold', position = position_stack(vjust = 0.5)) +
    scale_color_manual(values = order) + 
    scale_fill_manual(values = order) + 
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    theme_bw() +
    guides(color = 'none', fill = 'none', label = 'none') +
    theme(
      aspect.ratio = 6,
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
    ggsave(save_as, width = 2.5, height = 11)
}

figures_exac_vep <- function(given_vep = 'data/ExAc/VEP.csv') {
  vep <- 
    read_csv(
      given_vep, 
      col_types = cols_only(
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
    ) %>%
    distinct()
  
  vep %>%
    filter(mutation_type %in% c('Missense', 'Silent', 'Nonsense')) %>%
    count(mutation_type) %>%
    ungroup() %>%
    mutate(
      percent = n / sum(n),
      mutation_type = ordered(mutation_type, mutation_type[order(percent, decreasing = TRUE)])
    ) %>%
    ggplot(aes(y = percent, x = 'ExAc', fill = mutation_type)) +
    geom_col(color = 'black', alpha = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(-0.01, 1.01), expand = FALSE) +
    labs(x = '', y = '% of ExAc variants') +
    theme_bw() +
    theme(aspect.ratio = 5, panel.border = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank()) +
    ggsave('figures/ExAc-VEP-mutation-types.pdf')
}
