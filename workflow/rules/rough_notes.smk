

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