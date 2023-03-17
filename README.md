DOCUMENTATION:
Tool Installation:
To run the Python Pipeline, the following tools must be installed: (1) BowTie2 (to create an index for the HCMV genome and map donor reads to the resulting index); (2) SPAdes (to generate an assembly of the four patient transcriptomes); (3) BLAST+ (to query the nucleotide database for all members of the Betaherpesvirinae subfamily); and (4) Python 3.11.1 (to run Python Pipeline). A secure connection to the server via the terminal is also needed. 

Necessary Components:
The following .py files must be downloaded in order to run the code:
(1) PipelineProject.py (main driver code); and (2) Accession_to_FASTA.py (used to retrieve patient samples from NCBI)
Run the PipelineProject.py code to produce the PipelineProject.log containing the requested output from running the pipeline with all input reads. Both .py files should be in the same directory on your computer for the code to be executed. Make sure your working directory is set to the directory in which your .py files are located. A Python output log will be outputted with the following information specified in the Functionality section from the code. 

Testing Code with Sample Data:
Download the folder containing the sample data (4 truncated transcriptomes from the 2 patient donors), and move it into the same directory as the PipelineProject.py and Accession_to_FASTA.py files. Again, make sure your working directory is set to the directory in which your .py file are located. 

Functionality/Output:
HCMV transcriptomes from two patient donors from 2- and 6- days post-infection will be compared to all strains of the Human herpesvirus 5 (aka, Human cytomegalovirus or HCMV), as well as other virus strains. 
The following information will be outputted in the PipelineProject.log file:
(1) For each transcriptome, the number of read pairs mapped to the HCMV genome before and after Bowtie2 filtering will be outputted.
(2) In the SPAdes assembly of the four transcriptomes, the number of contigs with a length exceeding 1000 base pairs will be outputted.
(3) In the SPAdes assembly of the four transcriptomes, the number of base pairs in the assembly will be outputted.
(4) Using the longest contig from the SPAdes assembly, the subject accession, percent identity, alignment length, start of alignment in query, end of alignment in query, start of alignment in subject, end of alignment in subject, bit score, E-value, and subject title will be outputted for the top 10 hits of members of the Betaherpesvirinae subfamily. 
