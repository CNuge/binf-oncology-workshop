"""
 Systematic quality evaluation of copy number variants (CNVs) using IGV


note for lynch
PMS2 has pseudogenes that affect the ability to call variants with short reads
Need long reads to get reliable calls for the back half of the gene.



******
Lynch syndroms is an ideal case for applying long read

MLH1 is prone to epigenetic silencing as a disease mechanism
(promoter hypermethylation)

PMS2  has pseudogenes that affect the ability to call variants with short reads
Need long reads to get reliable calls for the back half of the gene.

2/5 most prone genes have problems solved by long read.
****




#### 
# Module 4


CNVs, in IGV - in whole exome looking at read depth for a region relative to nearby
- wouldn't apply in WGS as you're likely to get more sporadic coverage, higher variance across the genome.

What I see:
- lower depth in exons 6 and 7
- sporadic reads in the introns before the 6th and after the 7th exon, sharp drop off pre 6
- some sort of sharp drop in depth in the 8th exon - drops from 100 ->60 (~half)


- find a deletion in short read exome data in IGV, based on changes in coverage histogram
   sharp drop, then half depth of the controls, while maintaining the sequence quality.

"""