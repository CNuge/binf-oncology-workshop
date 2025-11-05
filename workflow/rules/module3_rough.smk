

ssh -i CBW.pem ubuntu@13.uhn-hpc.ca
conda activate dev
cd /home/ubuntu/CourseData/Module3

###########
# 
# victor martinez - iwk clinical genomics
# testing of heriditary breast and ovarian cancer
# about 10 genome analysts now; planning to add to team
# he attended/presented? at the conference tom was at

"""

For PGx calls they use the DRAGEN workflow to get start calls for the region.
DRAGEN has some reprocessing of these regions under the hood.

DRAGEN can check ploidy - this may break down in tumours


the tumour/normal is not applied for hematologic malignancies
    - cells are circulating in the blood/ purity >80% but you got to filter the germline cells


Free version of Franklin:
https://franklin.genoox.com/


For this workshop: Ensure “hg19” (GRCh37) is selected
(***shudders***)


Here's the RAD51c variant

https://franklin.genoox.com/clinical-db/variant/snp/chr17-58696863-A-G-hg38
from this one:
https://pubmed.ncbi.nlm.nih.gov/31782267/
I met Dr. Dawson on Tues evening.


They use phenotips to capture phenotype information.
- clinician notes are input to the patient metadata, then
- it gives the HPO terms for the phenotype


- people after free tools - https://varsome.com/ recommended by Amy, its becoming more
restricted with time with paywalls introduced

- note passing of the non annotated VCFs to the website; could likely do either.


- once vcfs loaded in, you can take the information and make your own classification
- links out to google schilar, mastermind, etc.

# free GPT one:
# https://chatgpt.com/g/g-niDu6xV5q-genomic-variant-interpretator?utm_source=gptshunter.com



- put in the genes with vcfs
- put in the hpo terms from phenotypes
- The “Phenotype Match” column will show scores indicating how well each variant’s gene matches the patient’s clinical features




"""


###########
#