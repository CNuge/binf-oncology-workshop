

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



#####
# Module 2 - intro to GATK


