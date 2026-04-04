import os
import ssl
from Bio import Entrez

ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = "2418821593@qq.com"

def fetch_sequence(gene_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
        record = handle.read()
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching sequence for {gene_id}: {e}")
        return None
    
def save_sequence_to_file(gene_id, content):
    try:
        os.makedirs("D:\\biocode\\3cytokine", exist_ok=True)
        filename = f"D:\\biocode\\3cytokine\\{gene_id}.fasta"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(content)
    except Exception as e:
        print(f"Error saving sequence to file {filename}: {e}")        
        
if __name__ == "__main__":
    gene_id = "NM_173923.2"
    sequence = fetch_sequence(gene_id)
    save_sequence_to_file(gene_id, sequence)
