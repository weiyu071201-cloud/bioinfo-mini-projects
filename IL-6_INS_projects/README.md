Project Overview

This project explores the variation of immune and metabolic genes across different species.
The immune gene is represented by IL-6, and the metabolic gene is insulin (INS).
By comparing CDS length and GC content, the project aims to investigate which type of gene is more conserved during evolution.


---

Workflow

grab.py: retrieves sequences from NCBI using accession numbers and saves them as FASTA files

analysis.py: extracts CDS regions, calculates CDS length and GC content, and stores results in a CSV file

graphic.py: visualizes the data from CSV and calculates standard deviation



---

Key Points

GC content is stored as a percentage string (e.g., “45.23%”), which needs to be converted to numeric format for plotting and analysis

When reading CSV files, pandas treats the first row as header by default, so extra care is needed if the file has no header



---

Problem & Solution

A key issue encountered was that ORF is not always equal to CDS.
In the mouse insulin 1 gene, there is an ORF located in the 5'-UTR region.
If simply selecting the first start codon (ATG), it leads to an incorrect 84 bp sequence.

To solve this, I referred to NCBI gene annotations and modified the algorithm to identify the longest ORF, which matches the annotated CDS.


---

Conclusion

In terms of CDS length, both IL-6 and INS genes show relatively small variation across species

In terms of GC content, the metabolic gene (INS) shows greater variation than the immune gene
