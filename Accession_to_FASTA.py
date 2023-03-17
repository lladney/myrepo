from Bio import Entrez              
Entrez.email = "lladney@luc.edu"

referenceSeq = (Entrez.efetch(db="nucleotide", rettype="fasta",
                retmode="text", id="NC_006273.2").read())
                                                                    
with open("referenceSeq_out.txt", "w") as outfile:        
    print(referenceSeq, file = outfile)
