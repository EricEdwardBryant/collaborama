read_lines("data/ClinVar/clinvar.vcf.gz", n_max = 100)

#Use read_vcf function from mutagenesis to read and write VCF tables
read_vcf("data/ClinVar/clinvar.vcf.gz", n_max = 100)
#In header ot VCF file, tells you reference file (GRCh38)
#Then reference column is reference genome

#str(.Last.value) #Look at structure of the last thing you looked at