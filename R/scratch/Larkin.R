

read_lines("data/ClinVar/clinvar.vcf.gz", n_max = 100)

read_vcf("data/ClinVar/clinvar.vcf.gz")

#first 26 lines are hashhas, dab deliminated file bunch of lines in the begining with hash marks. 
#to deal with that eric wrote a function called read_vcf. From the mutagenisis package. 
#mutagenisis::read_vcf

#this output, if we just look at the structure str(.last.value), use str to look at structure and 
# last.value is a variable that always exsist.

# contains chrom, pos, ID, ref, alt, qual. filter, (id is a clinvar ID, specific to clinvar). Filter = Pass to mark variants 
# passed a quality. filter and qual dont always get used.

# reference= in the intro to the vcf giveds you want reference genome.

data <- read_csv("data/Coriell/Catalog-Export.csv")

clear.nonsense <-
  data %>%
  filter(str_detect(Mutations, regex('[0-9][0-9]*?[*]|stop', ignore_case = T)),
  !str_detect(Mutations,"fs")
)

ID.clear <- clear.nonsense$ID


C.TO.T <-
  data %>%
  filter(str_detect(Mutations, regex('C>T$', ignore_case = T))
         
  )

ID.C.TO.T <- data$ID

write.csv(C.TO.T, "coriel-CtoT.csv")


CDS <- read.csv("data/CDS/Hsapiens-UCSC-hg38-validated.csv")

hg38 <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens

# check the protien for the noted mutation
sequence <-
CDS %>%
  filter(gene %in% c('TREX1', 'NPC1')) %>%
  arrange(tx, exon) %>%
  mutate(
    sequence = get_genomic_sequence(chr, strand, start, end, genome=hg38)
  ) %>%
group_by(tx, gene) %>%
  summarise(
    cds = str_c(sequence, collapse=""),
    protein = translate(cds)
  )


#position given in mutation  name times 3 is the start of the codon.

# position*3 to get the last, -2 for the first

sequence %>%
     filter(gene =="NPC1") %>%
    mutate(
      codon = str_sub(),
      
    )
  
  
  
  
  
  
  #then you can look at the actual codon, compare to the stop codons and determine what changes may have occured to make it a stop
  # por ejemplo they TCA could have become TGA or TAA to become a stome codon. The lab can only do the dark black arrows from the graphics
  # only when stop codon was the result of a c to T will they be to do perfect reversal.Still interested in nonesense to missense
  
  # MUST be near and NGG pam in order to edit. (or an NCC) has to be 10 to 16 (non inclusive) between the codon and the pam
# in order to be targetable.

