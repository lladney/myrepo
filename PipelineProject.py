# COMP483 Pipeline Project, Track 2 (Genome Assembly)
### STEP 1
outputLog = open("PipelineProject.log","w")                       # open and write to log file

import os                                                         # import os module
os.system("cd /home/lladney/COMP483")                             # change to working directory (COMP483 folder)
os.system("mkdir PipelineProject")                                # make a new directory for Pipeline Project
os.system("cd /home/lladney/COMP483/PipelineProject_Lara_Ladney") # change working directory to newly made directory

# retrieve transcriptomes from two patient donors from SRA
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045")

# convert SRA to paired-end FASTQ files
os.system("fastq-dump -I --split-files SRR5660030")
os.system("fastq-dump -I --split-files SRR5660033")
os.system("fastq-dump -I --split-files SRR5660044")
os.system("fastq-dump -I --split-files SRR5660045")

### STEP 2 (TRACK 2)
exec(open("Accession_to_FASTA.py").read())                        # retrieving the sequence from NCBI using Python code

os.system("bowtie2-build referenceSeq_out.txt  NC_006273.2")      # create an index for HCMV 

# save only reads that map to the HCMV index for use in assembly
os.system("bowtie2 -x NC_006273.2 -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S SRR5660030map.sam --al-conc SRR5660030_mapped_%.fq")

os.system("bowtie2 -x NC_006273.2 -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S SRR5660033map.sam --al-conc SRR5660033_mapped_%.fq")

os.system("bowtie2 -x NC_006273.2 -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S SRR5660044map.sam --al-conc SRR5660044_mapped_%.fq")

os.system("bowtie2 -x NC_006273.2 -1 SRR5660044_1.fastq -2 SRR5660045_2.fastq -S SRR5660045map.sam --al-conc SRR5660045_mapped_%.fq")

# write to log number of read pairs before and after mapping
# GET LINE COUNT BEFORE BOWTIE2 MAPPING
# Donor 1, 2 dpi
donor12Before = os.system("wc -l SRR5660030_1.fastq")             # get number of lines before filtering and assign to donor12Before variable
# os.system("wc -l SRR5660030_2.fastq")                           # CHECK: number of lines in raw fastq file 2 matches raw fastq file 1 (it does)
donor12Before_readpairs = (int(donor12Before)/4)                  # get number of read pairs by dividing number of lines by 4 and assign to donor12Before_readpairs variable 

donor12After = os.system("wc -l SRR5660030_mapped_1.fq")          # get number of lines before filtering and assign to donor12After variable
# os.system("wc -l SRR5660030_mapped_2.fq")                       # CHECK: number of lines in mapped fastq file 2 matches mapped fastq file 1 (it does)
donor12After_readpairs = (int(donor12After)/4)                    # get number of read pairs by dividing number of lines by 4 and assign to donor12After_readpairs variable

# Donor 1, 6dpi                                                                 
donor16Before = os.system("wc -l SRR5660033_1.fastq")             # get number of lines before filtering and assign to donor16Before variable
donor16Before_readpairs = (int(donor16Before)/4)                  # get number of read pairs by dividing number of lines by 4 and assign to donor16Before_readpairs variable 
donor16After = os.system("wc -l SRR5660033_mapped_1.fq")          # get number of lines before filtering and assign to donor16After variable
donor16After_readpairs = (int(donor16After)/4)                    # get number of read pairs by dividing number of lines by 4 and assign to donor16After_readpairs variable

# Donor 3, 2 dpi
donor32Before = os.system("wc -l SRR5660044_1.fastq")             # get number of lines before filtering and assign to donor32Before variable
donor32Before_readpairs = (int(donor32Before)/4)                  # get number of read pairs by dividing number of lines by 4 and assign to donor32Before_readpairs variable 
donor32After = os.system("wc -l SRR5660044_mapped_1.fq")          # get number of lines before filtering and assign to donor32After variable
donor32After_readpairs = (int(donor32After)/4)                    # get number of read pairs by dividing number of lines by 4 and assign to donor32After_readpairs variable

# Donor 4, 6 dpi
donor36Before = os.system("wc -l SRR5660045_1.fastq")             # get number of lines before filtering and assign to donor36Before variable
donor36Before_readpairs = (int(donor36Before)/4)                  # get number of read pairs by dividing number of lines by 4 and assign to donor36Before_readpairs variable 
donor36After = os.system("wc -l SRR5660045_mapped_1.fq")          # get number of lines before filtering and assign to donor36After variable
donor36After_readpairs = (int(donor36After)/4)                    # get number of read pairs by dividing number of lines by 4 and assign to donor36After_readpairs variable

# write print statements on # of read pairs before and after filtering for each donor/dpi to the log file
outputLog.write("Donor 1 (2dpi) had " + str(donor12Before_readpairs) + " read pairs before Bowtie2 filtering and " + str(donor12After_readpairs) + " read pairs after.")
outputLog.write("Donor 1 (6dpi) had " + str(donor16Before_readpairs) + " read pairs before Bowtie2 filtering and " + str(donor16After_readpairs) + " read pairs after.")
outputLog.write("Donor 3 (2dpi) had " + str(donor32Before_readpairs) + " read pairs before Bowtie2 filtering and " + str(donor32After_readpairs) + " read pairs after.")
outputLog.write("Donor 3 (6dpi) had " + str(donor36Before_readpairs) + " read pairs before Bowtie2 filtering and " + str(donor36After_readpairs) + " read pairs after.")

### STEP 3 (TRACK 2)
# assemble 4 transcriptomes together using SPAdes
os.system("spades.py -k 77,99,127 -t 4 --only-assembler --pe-1 1 SRR5660030_mapped_1.fq --pe-2 1 SRR5660030_mapped_2.fq --pe-1 2 SRR5660033_mapped_1.fq --pe-2 2 SRR5660033_mapped_2.fq --pe-1 3 SRR5660044_mapped_1.fq --pe-2 3 SRR5660044_mapped_2.fq --pe-1 4 SRR5660045_mapped_1.fq --pe-2 4 SRR5660045_mapped_2.fq -o HCMV_transcriptome_assembly/")

### STEP 4 (TRACK 2)
# change working directory to folder generated during SPAdes that contains contigs FASTA
os.system("cd ~/COMP483/PipelineProject/HCMV_transcriptome_assembly")

# write Python code to calculate the number of contigs with a length greater than 1000 bp in the assembly...
# and the number of bp in all of the contigs exceeding 1000 bp
from Bio.Seq import Seq                                          # imports sequence class
from Bio import SeqIO                                            # importing SeqIO class from BioPython to parse FASTA files
import numpy                                                     # import numpy module to get cumulative sums of contig lengths

handle = open("contigs.fasta")                                   # opening fasta file containing sequences

contigLenList = []                                               # create list to store contig lengths
for record in SeqIO.parse(handle, format = "fasta"):             # using a for loop to parse FASTA file 
                                                                 # for loop iterates through each contig
    recordLen = len(record)                                      # variable "recordLen" created to store contig lengths
    if recordLen > 1000:                                         # for loop iterates through each contig length and adds 1 to the
        contigLenList.append(recordLen)                          # append recordLen if length exceeds 1000 bp to countLenList
contigCount = len(contigLenList)                                 # take the length of the contigLenList to determine # of contigs exceeding 1000 bp
bpInContigs = sum(contigLenList)                                 # take the sum of all contig lengths stored in the list to determine # of bp in contigs
outputLog.write("There are " + str(contigCount) + " contigs > 1000 bp in the assembly.")
outputLog.write("There are " + str(bpInContigs) + " bp in the assembly.")

### STEP 5 (TRACK 2)
# write Python code to retrieve the longest contig from your SPAdes assembly  
maxRecord = 0                                                    # initialize maxRecord variable to 0
longestContig = ""                                               # create longestContig variable to store longest contig
for record in SeqIO.parse(handle, format = "fasta"):             # using a for loop to parse FASTA file 
                                                                 # for loop iterates through each contig
    recordLen = len(record)                                      # variable "recordLen" created to store contig lengths
    if recordLen > maxRecord:                                    # for loop iterates through each contig length and makes record
        recordLen = maxRecord                                    # append recordLen if length exceeds 1000 bp to countLenList
        longestContig = record.seq                               # make the sequence of the corresponding record the longestContig

# get all record sequences for Betaherpesvirinae virus
Entrez.email = "lladney@luc.edu"                                 # tell NCBI who you are
betaRecord = Entrez.esearch(db = "nucleotide", term = "Betaherpesvirinae")
betas = Entrez.read(betaRecord)
id_list = betas["IdList"]
id_string = ','.join(id_list)
betaSeqs = Entrez.efetch(db = "nucleotide", id = id_string, rettype = "fasta", retmode = "text")
output_fas = open("Betaherpesvirinae.fasta")
output_fas.write(betaSeqs)

# make a local datatbase of just sequences from the Betaherpesvirinae subfamily
os.system("makeblastdb -in Beta.fasta -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl")

# use the longest contig as BLAST+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily
outputFile = "betaherpesvirinaeResults.csv"                      # create outputFile to store results from blastn run

# create blastCommand to query just members of the Betaherpesvirinae subfamily
blastCommand = "blastn -query " + longestContig + "-db -Betaherpesvirinae -out " + outputFile + " -outfmt 6 sacc pident length qstart qend sstart send bitscore evalue stitle -max_hsps 1"
os.system("blastCommand")                                        # call the blastn command

# write top 10 hits to the log file (outputLog)
topTen = "head -n 10 outputFile | -out" + outputLog              # output first 10 lines of outputFile and send to outputLog
os.system("topTen")                                              # call the topTen command






