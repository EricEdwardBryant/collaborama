

read_lines("data/ClinVar/clinvar.vcf.gz", n_max = 100)

#first 26 lines are hashhas, dab deliminated file bunch of lines in the begining with hash marks. 
#to deal with that eric wrote a function called read_vcf. From the mutagenisis package. 
#mutagenisis::read_vcf

#this output, if we just look at the structure str(.last.value), use str to look at structure and 
# last.value is a variable that always exsist.

# contains chrom, pos, ID, ref, alt, qual. filter, (id is a clinvar ID, specific to clinvar). Filter = Pass to mark variants 
# passed a quality. filter and qual dont always get used.

# reference= in the intro to the vcf giveds you want reference genome.




