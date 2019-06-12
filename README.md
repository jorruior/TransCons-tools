## README

This readme is a guideline for any user that wants to use the main Methods in Ruiz-Orera et al: https://doi.org/10.1101/348326  


DEPENDENCIES: 

-**BLAST** (v2.3.0)

-**bedtools** (v2.28.0)

-**PAML**

-**prank** (v4)

-**Rfoot** (included in utils)

-**Python packages: Biopython, scipy**



1) **getRegions.py**: Extract the conserved and non-conserved regions for every gene in a GTF+FASTA, masking pseudogenes and addying a list of mouse-human orthologs from Ensembl Compara.
```
python3 getRegions.py --gtf <TRANSCRIPT GTF> --fasta <TRANSCRIPT FASTA> --out <OUT NAME> --db <BLAST DB> -O input/mmu_to_hsa_orthologs.txt -m input/mmu_pseudogenes.gtf
```


2) **getRNP.py**: Get a list of RNPs and coverage for a list of transcripts or regions:
```
python3 getRNP.py --input <TRANSCRIPT/REGIONS PRED> --sam <SAMFILE_PLUS,SAMFILE_REV/SAMFILE> --out <OUT NAME> --cds <BED CDS>
```


3) **featureCov.py**: Compute the overlap in the previously computed regions for a specific feature (BED file): RNA-seq and Ribo-seq, promoter, ORF, RNP, or CLIP-seq overlap.
```
python3 featureCov.py --input <REGIONS OUTPUT BED> -f <BED FEATURES> -o <OUT NAME> --stranded <yes/no>
```

**Guidelines for reproducibility of methods in Ruiz-Orera et al.:**

- The initial mouse transcript dataset corresponds to the annotated version in Ensembl v.89. Repeats were masked using RepeatMasker.

- The human transcript dataset used for building the database in BLAST was obtained from: https://figshare.com/articles/Ruiz-Orera_et_al_2017_/4702375

- The final list in Ruiz-Orera et al. was curated by eliminating lncRNA regions that had protein-coding orthologs in mouse (possible unannotated pseudogenes), putative misannotated UTR regions (located within 4kb from a sense protein-coding gene and/or with evidence of being part of the same gene using RNA-Seq data), or regions with a RNA-seq coverage lower than X.

- The full coordinates of translated ORFs can be found in the input folder (mmu89_t_orfs.fa/bed/gtf). getDNDS.py allows to compute dn/ds in a list of ORFs. Both species1 and species2 FASTA are needed (genomic alignment, two fasta files). For the article, the genomic alignments between mouse and human ORFs were used. The alignment of the 9 peptide candidates in lncRNAs can be reproduced:
```
python3 getDNDS.py -1 input/candidate_peptides_sp1.fa -2 input/candidate_peptides_sp2.fa -o candidates
```

- The folder 'tables' contain all raw data to reproduce the figures in the paper:
*Table 1 contains data for all considered regions (for Figs 1A,2B,3C,4B)*
*Table 2 contains data for equally sized gene subregions (for Fig 1B)*
*Table 3 contains Ribo-Seq coverage data for regions divided by read length (for Fig 3A)*
*Table 4 contains data for all considered genes (for Figs 2A,3B)*
