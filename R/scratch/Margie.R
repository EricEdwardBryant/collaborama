######################################
# 1. Group types of mutations and filter for genes that have >5 nonsense mutations
# 2. Classify genes based on COSMIC census (TSG, Oncogene, other)
# 3. Identify top genes with high proportion of nonsense mutations compared to the genome-wide proportion of (nonsense mutations)/(all mutations). 
#

library(dplyr)
library(tidyr)
library(data.table)

### Looking at COSMIC 

vep <- fread("data/COSMIC/VEP.csv")

#Summary of overall mutations
t_all <- length(vep$mutation_type) #12749594
m_all <- length(subset(vep$mutation_type, vep$mutation_type == "Missense")) #596436 #4.7% of mutations
nsm_all <- length(subset(vep$mutation_type, vep$mutation_type == "Nonsense")) #8444279 #66.2% of mutations
fs_all <- length(subset(vep$mutation_type, vep$mutation_type == "Frameshift")) #418217 #3.2%
s_all <- length(subset(vep$mutation_type, vep$mutation_type == "Silent")) #2991576 #23.4%
sp_all <- length(subset(vep$mutation_type, vep$mutation_type == "Splicing")) #298986 #2.34%

#Margie Make summary table of above

#First find subset of genes that have nonsense mutations
length(unique(vep$gene)) #18963 unique genes
genename <- c(unique(subset(vep$gene, vep$mutation_type == "Nonsense"))) #17508 genes have nonsense mutations

#subset to make loop faster
features = c("gene", "mutation_type")
vepselect = select(vep, features)

#save(vepselect, file = "data/COSMIC/VEP_subset_atleast5.RData")
#load("/Users/TinyDragon/github/collaborama/data/COSMIC/VEP_subset_atleast5.RData")

#Read file containing COSMIC Census
cosmic_genes <- read.csv("/Users/TinyDragon/github/collaborama/data/COSMIC/Census_allMon Nov 20 16_36_34 2017.csv")

##This data contains two "Tiers" of genes. Tier 1 contains established driver genes 
#and Tier 2 genes have less but emerging evidence supporting their role in cancer
#I elected to combine gene classes as below, and classify ambiguous and unclassified ones as "Other"

#Summary of cancer role status for each gene
table(cosmic_genes$Role.in.Cancer)
cosmic_genes2 <- cosmic_genes %>% select(Gene.Symbol, Role.in.Cancer, Molecular.Genetics) %>%
  filter(!Role.in.Cancer %in% c("", "fusion")) %>%
  rename("gene" = "Gene.Symbol") %>%
  mutate(Role.in.Cancer = as.vector(Role.in.Cancer)) %>%
  mutate(Role.in.Cancer = replace(Role.in.Cancer, Role.in.Cancer == "TSG, fusion", "TSG")) %>%
  mutate(Role.in.Cancer = replace(Role.in.Cancer, Role.in.Cancer == "oncogene, fusion", "oncogene")) %>%
  mutate(Role.in.Cancer = replace(Role.in.Cancer, Role.in.Cancer == "oncogene, TSG, fusion", "Other")) %>% #Not sure what to do with this category, named it "Other" for now, since it wouldnt
  mutate(Role.in.Cancer = replace(Role.in.Cancer, Role.in.Cancer == "oncogene, TSG", "Other")) 
table(cosmic_genes2$Role.in.Cancer)
#Code that can be used to merge the TSG,fusion and TSG categories, and oncogene,fusion and oncogene categories - should we do this?? or just name them all "other"?


#Counting and Grouping cosmic data by gene and mutation type
vepjoin <- left_join(vepselect, cosmic_genes2) %>%
  mutate(Role.in.Cancer = replace(Role.in.Cancer, is.na(Role.in.Cancer), "Other")) 

vepcircos <- table(Role = vepjoin$Role.in.Cancer, Mutation = vepjoin$mutation_type)
write.table(vepcircos, file = "/Users/TinyDragon/github/collaborama/data/COSMIC/vepcircos.txt")
#Use above table at http://mkweb.bcgsc.ca/tableviewer/visualize/ to make Circos visualization
#Non-contingency type ways to make these - for use in analysis. Above format necessary for Circos
#vepcircos <- count(vepjoin, Role.in.Cancer, mutation_type)
#vepcircos <- aggregate(Role.in.Cancer ~ mutation_type, data = vepjoin, FUN = sum)

vepsum <- count(vepselect, gene, mutation_type) 
# <- left_join(vepsum, cosmic_genes2) %>% #Commented out because I want to join it later
  #mutate(Role.in.Cancer = replace(Role.in.Cancer, is.na(Role.in.Cancer), "Other")) 

#Now transforming data to make types of mutations their own variables, gene by gene
vepsumwide  = gather(vepsum, variable, value, mutation_type) %>%
  unite(var, variable) %>% 
  group_by(gene) %>% 
  #mutate(group_row = 1:n()) %>%
  spread(value, n, fill = 0) %>% #having it fill the missings with 0 instead of NA since it just means the gene didn't have any of whatever type of mutation
  mutate(total_mut = Frameshift + Missense + Nonsense + Silent + Splicing)

#Finally, classify genes based on information from COSMIC
vepsumwide <- left_join(vepsumwide, cosmic_genes2) %>%
  mutate(Role.in.Cancer = replace(Role.in.Cancer, is.na(Role.in.Cancer), "Other")) 

vepsumwide <- mutate(vepsumwide, propnsm = Nonsense/total_mut) #make variable for proportion of nonsense to total mutation

#Test proportion of nonsense
#BUT how to test since one gene will be contributing the genome-wide rate? 
#Not really a two-sample test, since one sample is part of the other. Hmmm

#Just look at top?
vepsubset <- subset(vepsumwide, vepsumwide$Nonsense > 5) #filter for those with >5 NSM
vepsubset <- subset(vepsubset, vepsubset$propnsm > 0.078) #filter for proportion greater than the population mean estimated by genome-wide rate


###################################################
# Old loop from before I found fast dplyr way >_< #

#Now we want to find genes that have at least 5 nsm

p = c()
or = c()
length = c()
hugo = c()
total = c()
missense = c()
nsm = c()
frameshift = c()
silent = c()
splice = c()
start.time = Sys.time()
for (i in 100:102){   #test with a few
  #for (i in 1:length(genename)){  #for full set : took 49 minutes
  tryCatch({
    sub <- subset(vep, vep$gene == genename[i])
    if (length(subset(sub, sub$mutation_type == "Nonsense")) >=5){
      t <- length(sub$mutation_type) #12749594
      m <- length(subset(sub$mutation_type, sub$mutation_type == "Missense")) 
      nsm <- length(subset(sub$mutation_type, sub$mutation_type == "Nonsense")) 
      s <- length(subset(sub$mutation_type, sub$mutation_type == "Silent")) 
      ###Initial length and pos if using protein length and Amino acid position
      #Get initial lengths and positions - using PROTEIN instead of gene length
      hugo <- c(hugo, genename[i])
      print(i)
    } #end first while curly
    else {
      print(i)}
  }, #Trycatch end curly
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
} #end loop
end.time = Sys.time()
end.time - start.time #Took 49 minutes


#Trying to improve computing time


hugo = c()
start.time = Sys.time()
#for (i in 100:102){   #test with a few
for (i in 1:length(genename)){  #for full set : took 49 minutes
  tryCatch({
    sub <- subset(vep, vep$gene == genename[i])
    if (length(sub[,grep("Nonsense", sub$mutation_type)]) >=5){
      hugo <- c(hugo, genename[i])
      print(i)
    } #end first while curly
    else {
      print(i)}
  }, #Trycatch end curly
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
} #end loop
end.time = Sys.time()
end.time - start.time 


p = c()
or = c()
length = c()
#hugo = c()
total = c()
missense = c()
nonsense = c()
frameshift = c()
silent = c()
splice = c()
start.time = Sys.time()
for (i in 100:102){   #test with a few
  #for (i in 1:length(hugo)){  #for full set : took 49 minutes
  tryCatch({
    sub <- subset(vepselect, vepselect$gene == hugo[i])
    #if (length(sub$mutation_type[grep("Nonsense", sub$mutation_type)]) >=5){
    t <- length(sub$mutation_type) #12749594
    m <- length(sub$mutation_type[grep("Missense", sub$mutation_type)]) 
    nsm <- length(sub$mutation_type[grep("Nonsense", sub$mutation_type)]) 
    s <- length(sub$mutation_type[grep("Silent", sub$mutation_type)]) 
    fs <- length(sub$mutation_type[grep("Frameshift", sub$mutation_type)])
    sp <- length(sub$mutation_type[grep("Splicing", sub$mutation_type)])
    total = c(total, t)
    missense = c(missense, m)
    nonsense = c(nonsense, nsm)
    frameshift = c(frameshift, fs)
    silent = c(silent, s)
    splice = c(splice, sp)
    #hugo <- c(hugo, genename[i])
    print(i)
    # } #end first while curly
    # else {
    #   print(i)}
  }, #Trycatch end curly
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
} #end loop
end.time = Sys.time() #3.34 hours
end.time - start.time #Took 49 minutes

save(hugo, file = "data/COSMIC/cosmic_nsm_atleast5_hugo.RData") #File previously made with genes with 5 nsm, but found faster way to do it later
load("/Users/TinyDragon/github/collaborama/data/COSMIC/cosmic_nsm_atleast5_hugo.RData")


#########################################
#Old code from first meeting together
read_lines("data/ClinVar/clinvar.vcf.gz", n_max = 100)

#Use read_vcf function from mutagenesis to read and write VCF tables
read_vcf("data/ClinVar/clinvar.vcf.gz", n_max = 100)
#In header ot VCF file, tells you reference file (GRCh38)
#Then reference column is reference genome

#str(.Last.value) #Look at structure of the last thing you looked at

#Let's mess things up for real

corielldat <- read_csv("data/Coriell/Catalog-Export.csv")
#First task: find nonsense
coriell.narm <- na.omit(corielldat)

#The "*" denotes stop codon
coriell.nsm <- coriell.narm[grepl("\\*$", coriell.narm$Mutations),]
set.endstar <- coriell.nsm$ID

#C>T mutations are also of interest
coriell.nsm2 <- coriell.narm[grepl("C>T$", coriell.narm$Mutations),]
set.CT <- coriell.nsm2$ID

CDS <- read.csv("data/CDS/Hsapiens-UCSC-hg38-validated.csv")
hg38 <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
#To get sequence data, you need to know where in the genome, what chorm on and strand, and have ref genome

sequences <- 
  CDS %>%
  filter(gene %in% c('TREX1', 'NPC1')) %>%
  arrange(tx, exon) %>% #to stitch in the right order!!
  mutate(
    sequence = get_genomic_sequence(chr, strand, start, end, hg38) #creates sequence for each exome 
  ) %>%
  group_by(tx, gene) %>%
  summarise(
    cds = str_c(sequence, collapse = ''),
    aacid = translate(cds)
  )

#Here we see that the 738 position of NPC1 is a serine! 

#PAM is always NGG ABE editor in their lab. But we're REALLY looking for an NCC 12-16 ahead.
#CnnnnnnnnnnnnNGG
#CnnnnnnnnnnnnnnnnNGG
#AnnnnnnnnnnnnNGG

#Spacing to use to edit 
