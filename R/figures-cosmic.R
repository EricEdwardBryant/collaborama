figures_cosmic_vep <- function(given_vep = 'data/COSMIC/VEP.csv') {
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
    count(mutation_type) %>%
    ungroup() %>%
    mutate(
      percent = n / sum(n),
      mutation_type = ordered(mutation_type, mutation_type[order(percent, decreasing = TRUE)])
    ) %>%
    ggplot(aes(y = percent, x = 'COSMIC', fill = mutation_type)) +
    geom_col(color = 'black', alpha = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(-0.01, 1.01), expand = FALSE) +
    labs(x = '', y = '% of coding region mutations in COSMIC') +
    theme_bw() +
    theme(aspect.ratio = 5, panel.border = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank()) +
    ggsave('figures/COSMIC-VEP-mutation-types.pdf')
  
  vep %>%
    filter(mutation_type == 'Nonsense', str_length(ref_cds) == 3) %>%
    mutate(ref = toupper(ref_cds)) %>%
    count(ref) %>%
    ungroup() %>%
    mutate(
      percent = n / sum(n),
      ref = ordered(ref, ref[order(percent, decreasing = TRUE)])
    ) %>%
    ggplot(aes(y = percent, x = 'COSMIC', fill = ref)) +
    geom_col(color = 'black', alpha = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(-0.01, 1.01), expand = FALSE) +
    labs(x = '', y = '') +
    theme_bw() +
    theme(aspect.ratio = 5, panel.border = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank()) +
    ggsave('figures/COSMIC-VEP-nonsense-from-codon.pdf')
  
  vep %>%
    filter(mutation_type == 'Nonsense', str_length(alt_cds) == 3) %>%
    mutate(alt = toupper(alt_cds)) %>%
    count(alt) %>%
    ungroup() %>%
    mutate(
      percent = n / sum(n),
      alt = ordered(alt, alt[order(percent, decreasing = TRUE)])
    ) %>%
    ggplot(aes(y = percent, x = 'COSMIC', fill = alt)) +
    geom_col(color = 'black', alpha = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(-0.01, 1.01), expand = FALSE) +
    labs(x = '', y = '') +
    theme_bw() +
    theme(aspect.ratio = 5, panel.border = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank()) +
    ggsave('figures/COSMIC-VEP-nonsense-to-codon.pdf')
  
  vep %>%
    filter(mutation_type == 'Nonsense', str_length(alt_cds) == 3, str_length(ref_cds) == 3) %>%
    mutate(change = paste(toupper(ref_cds), 'to', toupper(alt_cds)) ) %>%
    count(change) %>%
    ungroup() %>%
    mutate(
      percent = n / sum(n),
      change = ordered(change, change[order(percent, decreasing = TRUE)])
    ) %>%
    ggplot(aes(y = percent, x = 'COSMIC', fill = change)) +
    geom_col(color = 'black', alpha = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(-0.01, 1.01), expand = FALSE) +
    labs(x = '', y = '') +
    theme_bw() +
    theme(aspect.ratio = 5, panel.border = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank()) +
    ggsave('figures/COSMIC-VEP-nonsense-from-codon-to-codon.pdf')
}

figures_cosmic_consequence <- function() {
  data <- 
    read_tsv('data/COSMIC/CosmicMutantExport.tsv.gz') %>%
    filter(
      GRCh == 38, 
      SNP == 'y', !is.na(SNP),
      str_detect(`Mutation Description`, 'Substitution'),
      !is.na(`FATHMM prediction`)
    ) %>%
    select(
      coord = `Mutation genome position`, 
      strand = `Mutation strand`, 
      gene = `Gene name`, 
      cosmic_id = `Mutation ID`, 
      CDS = `Mutation CDS`, 
      AA = `Mutation AA`, 
      primary_site = `Primary site`,
      description = `Mutation Description`, 
      FATHMM = `FATHMM prediction`
    ) %>%
    mutate(
      description = case_when(
        str_detect(description, 'Missense') ~ 'Missense',
        str_detect(description, 'Nonsense') ~ 'Nonsense',
        str_detect(description, 'silent')   ~ 'Silent'
      ),
      FATHMM = case_when(
        FATHMM == 'NEUTRAL' ~ 'Neutral',
        FATHMM == 'PATHOGENIC' ~ 'Pathogenic'
      )
    )
  
  summary <- 
    count(data, description, FATHMM) %>%
    group_by(description) %>%
    mutate(
      total = sum(n),
      percent = n / total
    ) %>%
    ungroup()
  
  summary %>%
    mutate(
      description = ordered(description, levels = c('Nonsense', 'Missense', 'Silent'))
    ) %>%
    ggplot(aes(x = description, y = percent, fill = FATHMM)) +
    geom_col(color = 'black', alpha = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c('grey', 'red')) +
    labs(x = '', y = 'COSMIC coding single nucleotide variants') +
    coord_flip(ylim = c(-0.01, 1.01), expand = FALSE) +
    theme_bw() +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(aspect.ratio = 0.5, panel.border = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.y = element_blank(), legend.direction = 'horizontal', legend.position = 'top') +
    ggsave('figures/COSMIC-FATHMM-SNV-consequences.pdf', height = 3, width = 6)
      
  summary <-
    count(data, description, primary_site) %>%
    group_by(primary_site) %>%
    mutate(
      total = sum(n),
      percent = n / total
    ) %>%
    filter(all(n > 100)) %>%
    ungroup()
  
  summary %>%
    mutate(
      description = ordered(description, levels = c('Silent', 'Missense', 'Nonsense')),
      primary_site = ordered(primary_site, levels = primary_site[description == 'Nonsense'][order(percent[description == 'Nonsense'])])
    ) %>%
    ggplot(aes(x = primary_site, y = percent, fill = description)) +
    geom_col(color = 'black', alpha = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c('grey', 'blue', 'red')) +
    labs(x = '', y = 'COSMIC coding single nucleotide variants') +
    coord_flip(ylim = c(-0.01, 1.01), expand = FALSE) +
    theme_bw() +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(aspect.ratio = 0.2, panel.border = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.y = element_blank(), legend.direction = 'horizontal', legend.position = 'top') +
    ggsave('figures/COSMIC-primary-site-SNV.pdf', width = 8, height = 4)
    
  summary %>%
    filter(description == 'Nonsense') %>%
    mutate(
      description = ordered(description, levels = c('Silent', 'Missense', 'Nonsense')),
      primary_site = ordered(primary_site, levels = primary_site[description == 'Nonsense'][order(percent[description == 'Nonsense'])])
    ) %>%
    ggplot(aes(x = primary_site, y = percent, fill = description)) +
    geom_col(color = 'black', alpha = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c('red')) +
    labs(x = '', y = 'COSMIC coding single nucleotide variants') +
    coord_flip(ylim = c(0, 0.05), expand = FALSE) +
    theme_bw() +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(aspect.ratio = 0.2, panel.border = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.y = element_blank(), legend.direction = 'horizontal', legend.position = 'top') +
    ggsave('figures/COSMIC-primary-site-SNV-nonsense-only.pdf', width = 8, height = 4)
}
