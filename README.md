# sc-rna-seq-mini-demo
This repository is a demo workflow for single-cell RNA analysis with R

## sc-RNA micro demo RNotebook
This notebook contains simplest tutorial with the PBMC 3k data set, process includes:
- QC and filtering
- normalization 
- scaling, PCA, and clustering 
- find marker genes
- then, since cluster 3 is very distant from other clusterings, evaluate DE genes in cluster 3

## llm sequence analysis of the gene regulatory regions
- use python code to fetch 1kb upstream from the DE genes selected in cluster 3 and get a fasta file
- Step-by-Step LLM Pipeline with Hugging Face and PyTorch
    -  Step 1: Load FASTA sequences
    -  Step 2: Tokenize using 6-mer (for DNABERT) or character-based (for Nucleotide Transformer)
    - Step 3: Run through pretrained model to get embeddings
    - Step 4: Optional classification — enhancer/promoter inference
        - download enhancer from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) 
| Option            | Setting                                     |
| ----------------- | ------------------------------------------- |
| **clade**         | Mammal                                      |
| **genome**        | Human                                       |
| **assembly**      | hg19 (FANTOM5 enhancer track is under hg19) |
| **group**         | Regulation                                  |
| **track**         | FANTOM NET Enhancers                           |
| **table**         | enhancers                                   |
| **region**        | genome                                      |
| **output format** | BED - browser extensible data               |
  
        - download TSS from [FANTOM](https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/), file:hg19.cage\_peak\_phase1and2combined\_coord.bed.gz 
        - extract top 50 lines from each files as subsets to create references set for the demo purpose
        - download reference genome hg19.fa.gz from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)
        - then, perform a zero-shot classification based on cosine similarity

## another use case   
