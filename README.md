# VCF-parsing
This is a collection of useful scripts for parsing VCF files and extracting data for downstream analyses.

**extract-GT-V2-1-exclude-miss-called-positions-give-GT-too.pl**. 

By this script you can extract genotypes in a vcf file in two formats:  
"fasta" format  
#>ind1-1.  
hap-1.  
#>ind1-2.  
hap-2.  
"tab" format makes a matrix for haplotypes:  
ind1-1	h	a	p	-	1.  
ind1-2	h	a	p	-	2.   


Here is how you run it:   
extract-GT-V1.pl -format fasta/tab -vcf SNP.vcf -phased T/F -missingGT T/F. 



