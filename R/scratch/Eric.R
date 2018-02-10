read_csv('data/Coriell/Catalog-Export.csv') %>%
  filter(str_detect(Mutations, regex('[*]|TER|C>T', ignore_case = TRUE))) %>%
  select(ID, Gene, Mutations) %>%
  write_csv('data/Coriell/coriell-nonsense.csv')

curated <- read_csv('data/Coriell/coriell-nonsense-curated.csv')

nonsense <-
  curated %>%
  mutate(gene = case_when(
    gene == 'CKN1'   ~ 'ERCC8',
    gene == 'ARH'    ~ 'LDLRAP1',
    gene == 'CLN2'   ~ 'TPP1',
    gene == 'RECQL2' ~ 'WRN', # maybe BLM?
    gene == 'RECQL3' ~ 'BLM', # maybe WRN?
    TRUE ~ gene
  )) %>%
  group_by(gene, aa, aa_coord, ct_coord) %>%
  summarise(IDs = str_c(ID, collapse = ' | ')) %>%
  ungroup() %>%
  select(IDs, everything()) %>%
  mutate(
    cds_coord = ifelse(!is.na(aa_coord), (aa_coord * 3L) - 2L, (ceiling(ct_coord / 3) * 3) - 2L),
    aa_coord  = ifelse(!is.na(aa_coord), aa_coord, ceiling(ct_coord / 3)),
    cds_start = cds_coord,
    cds_end   = cds_coord
  )

CDS <- 
  read_csv('data/CDS/Hsapiens-UCSC-hg38-validated.csv') %>%
  filter(gene %in% nonsense$gene) %>%
  add_exon_details() %>%
  mutate(
    cds_start = exon_position,
    cds_end   = (exon_position + exon_width) - 1L
  )

joined <- 
  fuzzyjoin::genome_left_join(nonsense, CDS, c('gene', 'cds_start', 'cds_end')) %>%
  mutate(
    exon_boundary_dist = cds_coord - exon_position,
    codon_start = ifelse(strand == '+', start + exon_boundary_dist,      end - (exon_boundary_dist + 2L)),
    codon_end   = ifelse(strand == '+', start + exon_boundary_dist + 2L, end - exon_boundary_dist),
    sequence = get_genomic_sequence(chr, strand, codon_start, codon_end, hg38),
    aa_check = translate(sequence)
  ) %>%
  filter(aa == aa_check | is.na(aa)) %>%
  mutate(
    sequence_context = str_c(
      get_genomic_sequence(chr, strand, codon_start - 30, codon_start - 1L, hg38),
      get_genomic_sequence(chr, strand, codon_start, codon_end, hg38) %>% tolower(),
      get_genomic_sequence(chr, strand, codon_end + 1L, codon_end + 30, hg38)
    )
  ) %>%
  select(
    IDs, tx, gene = gene.x, aa = aa_check, ct_coord, aa_coord, cds_coord, 
    chr, strand, codon_start, codon_end, codon = sequence, sequence_context
  ) %>%
  group_by_at(vars(-tx)) %>%
  summarise(tx = str_c(tx, collapse = ' | ')) %>%
  select(IDs, tx, everything())

write_csv(joined, 'data/Coriell/coriell-nonsense-with-sequence-context.csv')
  
hg38 <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens

sequences <-
  CDS %>%
  filter(gene %in% c('TREX1', 'NPC1')) %>%
  arrange(tx, exon) %>%
  mutate(
    sequence = get_genomic_sequence(chr, strand, start, end, hg38)
  ) %>%
  group_by(tx, gene) %>%
  summarise(
    cds = str_c(sequence, collapse = ''),
    aa  = translate(cds)
  )

sequences %>%
  filter(gene == 'NPC1') %>%
  mutate(
    codon = str_sub(cds, 2212, 2214),
    upstream = str_sub(cds, 2214 - 20, 2214)
  ) %>%
  View()
