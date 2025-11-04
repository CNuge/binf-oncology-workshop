

######
# server login
ssh -i CBW.pem ubuntu@13.uhn-hpc.ca

#installed my own dev env
conda activate dev


(base) ubuntu@ip-10-0-1-243:~/CourseData$ ls
Module2  Module3  Module4  tools

(base) ubuntu@ip-10-0-1-243:~$ ls workspace
test.txt  test2.txt

# I can write files in workspace?

####
# Docker info

which docker


docker image list


##############################
# Module 2
### whats the data?
(dev) ubuntu@ip-10-0-1-243:~/CourseData/Module2/BAMs$ ls
C-GIAB.normal.regions_of_interest.bam
C-GIAB.normal.regions_of_interest.bam.bai
C-GIAB.tumour.regions_of_interest.bam
C-GIAB.tumour.regions_of_interest.bam.bai
(dev) ubuntu@ip-10-0-1-243:~/CourseData/Module2/BAMs$ pwd
/home/ubuntu/CourseData/Module2/BAMs

giab tumour/normal paired data

#intersting... I guess we can call variants by comparing the
#tumor to normal to see what has changed - neat variation.


# Dr. Gaston - work at IWK, clinical genomics
# did PhD in rare disease research.


###
# intrepretation of cancer genomics

# - Intrepretation of cancer-associated mutations is always changing
# more dynamic due to the evolution of new therapies. Things that used to indicate
# big problems may now indicate that a tumour is druggable

# clinical tests need to be certified in order to get into the clinic
# 


# tumour sequencing - 70-100x (more variance, harder to capture?)
# says they use PCR - really hope there is some pre-align cleaning that wasn't mentioned.


# Questions:
# joint calling within a cancer setting? (germline v cohort, possibly with a tumour sample)
# Long read not mentioned - would you use it with infinitie money?
# Long read - how important are structural variants within cancer characterization and treatment?


# Cancer data
# not just 0/0 0/1 1/1
# lots of different mutations may be present
# a tumour is a mix of normal or tumour cells.
# subclonal mutations and tumour heterogeneity
# 
# ^this is why you'd have higher depth in the tumor sequencing
# basically trying to sample a population.

# tumour - normal matched clinical sequencing.


# panel of normals - instead of sequencing every normal, take 40-60 normal people as a reference cohort
# can pass this into mutec2 instead of a matched normal.

# Iwk - use nextflow - older stuff was built in toil: https://toil.readthedocs.io/en/latest/
# whenever you modify a test / binf workflow it needs to be clinically recertified.
# therefore they are a little slower to update stuff.

# for variant intrepreatation - they use Illumina imagen
# I asked if intrepreter had experience with Golden Helix and such
# they said they were fine - 
# con with the illumina tool is you're pretty locked in on the front end.

# Mutect2? - in gatk
# https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
# call variants in a tumour-normal matched fashion


#####
# Module 2 - intro to GATK

######
# Docker

#Q1. Question: What is the path for the docker executable installed on your system?
 which docker
#/usr/bin/docker

#Q2. 
#Question: What are the names of the 3 docker images on the system? What repository did they come from? Can you use this informaton to find where on the internet the repository is for the two images that come from the same repository?
    docker image list
#REPOSITORY                                        TAG                   IMAGE ID       CREATED       SIZE
#public.ecr.aws/aws-genomics/broadinstitute/gatk   4.2.6.1-corretto-11   9c0cb40e920b   2 years ago   5.01GB
#pegi3s/fastqc                                     latest                9482f5c46ba2   2 years ago   579MB
#pegi3s/bwa                                        latest                d970a6e5eeaf   5 years ago   96.6MB

# can run the command as such:     docker run pegi3s/bwa bwa

#Q3 Question: What is the version of BWA being run with this docker container?

#(dev) ubuntu@ip-10-0-1-243:~/CourseData/Module2/FastQs$ docker run pegi3s/bwa bwa
#Program: bwa (alignment via Burrows-Wheeler transformation)
#Version: 0.7.17-r1188

#####
# me irl:
conda config --add channels bioconda
conda install gatk4

cd /home/ubuntu/CourseData/Module2/FastQs

zless normal_test_1.fastq.gz
zcat normal_test_1.fastq.gz
# old timey fastqs with the repeat for the header line!
zless normal_test_2.fastq.gz

# Q4. Questions: How many base pairs are the sequencing reads? What is the instrument ID of the sequencer these sequencing reads were generated on? What sequencing lane were the reads from in the normal and the tumour? What were the flowcell IDs?
zcat normal_test_1.fastq.gz | wc -l | awk '{print $1 / 4}'
#100K reads

zcat normal_test_1.fastq.gz | head
zcat tumour_test_1.fastq.gz | head
# What is the instrument ID of the sequencer these sequencing reads were generated on? 
# SRR30646153.1 I think?

#What sequencing lane were the reads from in the normal and the tumour? 
# What were the flowcell IDs?
# normal
# H73FJDSXC - but just cuz I know that looks like a lane id
# tumor
# HK25HDSX7 - different flowcell

#what lane did it come from? 
# normal
# lane 4
# tumor
# lane 3

#400000 / 4
zgrep "@" normal_test_1.fastq.gz | wc -l
#Question:* Why would this be a bad idea in practice?
# - @ is a valid phred encoding - so you would count those lines in a lower qual file

zcat normal_test_2.fastq.gz | wc -l
#400000 / 4

zcat normal_test_1.fastq.gz  | awk 'NR%4==2 {print length($0)}' > all_lens.txt

# Q5. Questions: How many reads are in the files? Are they all the same size?
#100k reads
# yes
ipython

    all_same = True
    with open("all_lens.txt") as file:
        for i, line in enumerate(file):
            line_len = line.rstrip()
            if int(line_len) != 150:
                all_same = False
                print(f"the {i} read is length: {line_len}")
    print("all the same?")
    print(all_same)


### 
# BED files
cd /home/ubuntu/CourseData/Module2/accessory_files
cat gene_target_regions.bed

## Q6. Question: What does each column represent?
# CHR chromStart chromEnd
# 



## Q7. Question: What combination of characters are used to represent tabs and newlines?
cp gene_target_regions.bed test.bed
nano test.bed
# \t and \n

cat gene_target_regions.bed
cat -A gene_target_regions.bed
# ^ and $ in the output of this
# prints it in regex


###
# BAM files

# cd /home/ubuntu/CourseData/Module2/BAMs
samtools head C-GIAB.normal.regions_of_interest.bam

samtools head C-GIAB.normal.regions_of_interest.bam | less
#ref with alts and decoys

# Q8. Question: What are the sort orders of the sample BAM files?
samtools head C-GIAB.normal.regions_of_interest.bam | grep "SO"
# coordinate sort order
# @HD     VN:1.4  SO:coordinate

# Q9. Question: What programs were used on this file?
samtools head C-GIAB.normal.regions_of_interest.bam | grep "SQ"

samtools head C-GIAB.normal.regions_of_interest.bam | grep "PG"

samtools head C-GIAB.normal.regions_of_interest.bam | tail
# last line ,dragen, samtools
# bunch on the previous line.
# @PG     ID:samtools     PN:samtools     PP: DRAGEN SW build     VN:1.12 CL:samtools view -b -@ 4 -M -L gene_target_regions.bed C-GIAB.bam


# Q10. Question: What do the -n, -@ and -O command-line flags do?

samtools -h
 -n         Sort by read name (natural): cannot be used with samtools index
 -@, --threads INT
 -O, --output-fmt 

# use samtools sort order so that it is sorted by read name instead of coordinate sort order
samtools sort -n -@ 4 -O BAM -o outtest.bam C-GIAB.normal.regions_of_interest.bam


# Q11 Question: If you use zless to look at these FastQ files, what do you notice compared to the ones you looked at previously?

#TODO - this fails when I pipe it:
samtools sort -n -@ 4 -O BAM -o outtest.bam C-GIAB.normal.regions_of_interest.bam | \
samtools view -h -u | \
samtools fastq -@ 2 -1 test_out.1.fastq.gz -2 test_out.2.fastq.gz -

#
samtools fastq -1 test_out_retry.1.fastq.gz -2 test_out_retry.2.fastq.gz -@ 2 outtest.bam
zless test_out_retry.1.fastq.gz
# difference is the lack of information on the header and the + row - we have lost some of the metadata

######
# vcfs

cd ~/CourseData/Module2/VCFs
ls
# Q12. Question: What file format specification was used for each file? Is it the same or different?
bcftools head mutect2_tumour_normal.filtered.regions_of_interest.vcf.gz | less
##fileformat=VCFv4.2
head test_panel.GRCh38.minimal.vcf
##fileformat=VCFv4.2

# Q13. Question: In the Mutec2 generated VCF how many different values can the FILTER field take?
bcftools head mutect2_tumour_normal.filtered.regions_of_interest.vcf.gz | grep "FILTER" | wc -l 
#23 -1 for the head ergo
#22 different listed

# Q14. Question: What is the ID in the FORMAT field for the Approximate read Depth at a site/variant?
bcftools head mutect2_tumour_normal.filtered.regions_of_interest.vcf.gz | grep "depth" 
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">

# Q15. Question: What are the genomic coordinates of the first variant listed in the file? 
#What are the Reference and Alt alleles? What entries are given in the filter field?
bcftools head mutect2_tumour_normal.filtered.regions_of_interest.vcf.gz -n 1
#chr9    21978178        .       T       G       .       base_qual;normal_artifact;strand_bias;weak_evidence       AS_FilterStatus=weak_evidence,base_qual,strand_bias;AS_SB_TABLE=31,36|0,5;DP=76;ECNT=1;GERMQ=93;MBQ=37,11;MFRL=402,447;MMQ=60,60;MPOS=48;NALOD=1.54;NLOD=10.44;POPAF=6;TLOD=3.44  GT:AD:AF:DP:F1R2:F2R1:FAD:SB      0/0:42,3:0.029:45:14,0:15,0:41,2:19,23,0,3      0/1:25,2:0.106:27:7,2:9,0:25,2:12,13,0,2

# filtered on: base_qual;normal_artifact;strand_bias;weak_evidence



## Misc notes

#cd /home/ubuntu/CourseData/Module2/BAMs
#samtools stats C-GIAB.tumour.regions_of_interest.bam

# fun fact
# bed file thickStart thickEnd itemRgb - used by IGV to draw lines on genome.

# Recall:
#    TTA[G/GC]CCCC
# - BAM, BED files are 0 based
#   half-open rangers. start in the range, end is not
#   3
# - SAM, VCF are 1-based 
#   Sam and BAM are slightly different
#   closed both start and end included in the range
#  4  G  GC


# data is a GiaB tumor sample.
# https://www.nature.com/articles/s41597-025-05438-2
# interesting broadly consented sample.






