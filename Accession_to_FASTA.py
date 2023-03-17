from Bio import Entrez                                          # import Entrez module to retrieve sequences from NCBI
Entrez.email = "lladney@luc.edu"                                # tell NCBI who you are

referenceSeq = (Entrez.efetch(db="nucleotide", rettype="fasta", # fetch reference sequence using efetch()
                retmode="text", id="NC_006273.2").read()) 
                                                                    
with open("referenceSeq_out.txt", "w") as outfile:              # open referenceSeq_out.txt as an outfile to write to
    print(referenceSeq, file = outfile)                         # print the reference sequence
