analyses$BRCA_summary <- function() {
  ClinVar_Model_CtoT <- read_csv("data/Guides/ClinVar-Model-CtoT.csv")
  # Note that there are 5 variants annotated as BRCA1:672|NBR2:10230. none are targetable by ABE or BE3
  ClinVar_Model_CtoT %>% filter(GENEINFO %in% c('BRCA1:672', 'BRCA2:675')) -> BRCA
  BRCA_stats <- table(transmute(BRCA, Gene = GENEINFO, BE3 = has(BE3), ABE = has(ABE), nonsense = str_detect(MC, 'nonsense')))
  reshape2::melt(BRCA_stats) %>% arrange(Gene, desc(nonsense), desc(BE3), desc(ABE))
  
  apply(BRCA_stats, 'nonsense', sum)
  apply(BRCA_stats, 'BE3',      sum)
  apply(BRCA_stats, 'ABE',      sum)
  margin.table(BRCA_stats, 1)
  margin.table(BRCA_stats, 2)
  margin.table(BRCA_stats, 3)
  
  table(ClinVar_Model_CtoT %>% transmute(BE3 = has(BE3), ABE = has(ABE)))
}