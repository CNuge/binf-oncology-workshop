
cd /home/ubuntu/CourseData/Module2

ls FastQs

#### 
# FastQC

#docker run --rm -v /home/ubuntu/CourseData:/data -v /home/ubuntu/workspace:/workspace pegi3s/fastqc -o /workspace/FastQC/ /data/Module2/FastQs/normal_test_1.fastq.gz


# 
#Question: How many reads are in the normal_test_1.fastq.gz file according to its FastQC output? How many total bases?
#100000
#Question: What is the GC % of th data?
#41
#Question: What is the average quality per read?
#36


#####
# FastQC

fastqc FastQs/normal_test_1.fastq.gz FastQs/normal_test_2.fastq.gz FastQs/tumour_test_1.fastq.gz FastQs/tumour_test_2.fastq.gz


lr FastQs/
cp FastQs/*.html ~/workspace


#####
# alignment
#Question: What other files are located in the ~/CourseData/tools/reference/GATK-bundle-GRCh38/ directory? Hint: The various iterations on Homo_sapiens_assembly38.fasta with extra extensions are index files used by various programs to efficiently use the FASTA reference file.
ls ~/CourseData/tools/reference/GATK-bundle-GRCh38/ 
# its the index files!

Create the command-line, using the docker elements and correct paths for docker to align the two normal FastQ files to the provided Human reference, and output to a SAM file in the ~/CourseData/Module2/Processing

/home/ubuntu/CourseData/Module2/FastQs/normal_test_1.fastq.gz
/home/ubuntu/CourseData/Module2/FastQs/normal_test_2.fastq.gz

/home/ubuntu/CourseData/Module2/FastQs/tumour_test_1.fastq.gz
/home/ubuntu/CourseData/Module2/FastQs/tumour_test_2.fastq.gz

bwa mem ~/CourseData/tools/reference/GATK-bundle-GRCh38/Homo_sapiens_assembly38.fasta  /home/ubuntu/CourseData/Module2/FastQs/tumour_test_1.fastq.gz /home/ubuntu/CourseData/Module2/FastQs/tumour_test_2.fastq.gz -o ~/CourseData/Module2/Processing/tumour_test_1.sam

bwa mem ~/CourseData/tools/reference/GATK-bundle-GRCh38/Homo_sapiens_assembly38.fasta  /home/ubuntu/CourseData/Module2/FastQs/normal_test_1.fastq.gz /home/ubuntu/CourseData/Module2/FastQs/normal_test_2.fastq.gz -o ~/CourseData/Module2/Processing/tumour_test_1.sam



samtools sort -@ 4 --write-index -O bam -o ~/CourseData/Module2/Processing/tumour_test_1.bam ~/CourseData/Module2/Processing/tumour_test_1.sam
samtools sort -@ 4 --write-index -O bam -o ~/CourseData/Module2/Processing/normal_test_1.bam ~/CourseData/Module2/Processing/normal_test_1.sam

#Question: What command-line flag could you give to the samtools sort command above to have it automatically create an index file for the output BAM?
--write-index 

# Question: If you wanted to sort the SAM, but output to a SAM file instead how would you modify the file?
# change the -O to a sam
# -u  -l 0
# What about if you wanted to create a BAM file but change its compression level to 0? 
# samtools view -h -b -l 0 input.bam -o output.bam
# Making a binary but uncompressed file?
# 


####
# gatk commands

docker run -v /home/ubuntu/CourseData:/gatk/CourseData --rm public.ecr.aws/aws-genomics/broadinstitute/gatk:4.2.6.1-corretto-11 gatk \
AddOrReplaceReadGroups -I /gatk/CourseData/Module2/Processing/normal.bam -O /gatk/CourseData/Module2/Processing/normal.rg.bam --RGID normal_lane1 --RGLB normal_lane1 --RGSM normal --RGPL ILLUMINA --RGPU Illumina

docker run -v /home/ubuntu/CourseData:/gatk/CourseData --rm public.ecr.aws/aws-genomics/broadinstitute/gatk:4.2.6.1-corretto-11 gatk \
SetNmMdAndUqTags -I /gatk/CourseData/Module2/Processing/normal.rg.bam -O /gatk/CourseData/Module2/Processing/normal.fixed.bam -R /gatk/CourseData/tools/reference/GATK-bundle-GRCh38/Homo_sapiens_assembly38.fasta

docker run -v /home/ubuntu/CourseData:/gatk/CourseData --rm public.ecr.aws/aws-genomics/broadinstitute/gatk:4.2.6.1-corretto-11 gatk \
MarkDuplicates -I /gatk/CourseData/Module2/Processing/normal.fixed.bam -O /gatk/CourseData/Module2/Processing/normal.dup_marked.bam -M /gatk/CourseData/Module2/Processing/normal.duplicate_metrics.txt

docker run -v /home/ubuntu/CourseData:/gatk/CourseData --rm public.ecr.aws/aws-genomics/broadinstitute/gatk:4.2.6.1-corretto-11 gatk \
BaseRecalibrator -I /gatk/CourseData/Module2/Processing/normal.dup_marked.bam -R /gatk/CourseData/tools/reference/GATK-bundle-GRCh38/Homo_sapiens_assembly38.fasta --known-sites /gatk/CourseData/tools/reference/GATK-bundle-GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /gatk/CourseData/tools/reference/GATK-bundle-GRCh38/af-only-gnomad.hg38.vcf.gz -O /gatk/CourseData/Module2/Processing/normal.bqsr_recal.table

docker run -v /home/ubuntu/CourseData:/gatk/CourseData --rm public.ecr.aws/aws-genomics/broadinstitute/gatk:4.2.6.1-corretto-11 gatk \
ApplyBQSR -R /gatk/CourseData/tools/reference/GATK-bundle-GRCh38/Homo_sapiens_assembly38.fasta -I /gatk/CourseData/Module2/Processing/normal.fixed.bam --bqsr-recal-file /gatk/CourseData/Module2/Processing/normal.bqsr_recal.table -O /gatk/CourseData/Module2/BAMs/normal.recalibrated.bam