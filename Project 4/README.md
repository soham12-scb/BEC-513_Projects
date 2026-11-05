# Project 4: GDC RNA-Seq Gene Expression Dataset

## Overview
**Project 4** contains a collection of RNA-Seq data files downloaded from the Genomic Data Commons (GDC). This dataset consists of gene expression counts processed using the STAR aligner and GENCODE annotations.

## Project Structure
The project is structured as a series of case-specific directories. Each directory typically contains:
- **`*.rna_seq.augmented_star_gene_counts.tsv`**: The primary gene expression count file.
- **`logs/`**: Directory containing metadata and processing logs (`*.parcel` files).
- **`annotations.txt` (Optional)**: Clinical or case-level annotations associated with the sample (e.g., prior malignancies, clinical status).

## Data Format
The expression files are tab-separated value (TSV) files (`augmented_star_gene_counts.tsv`) with the following structure:

- **Header**: Contains metadata such as the gene model version (e.g., GENCODE v36).
- **Columns**:
    - `gene_id`: Ensembl gene identifier.
    - `gene_name`: Official gene symbol.
    - `gene_type`: Type of gene (e.g., protein_coding).
    - `unstranded`, `stranded_first`, `stranded_second`: Raw read counts using different strand-specific library strategies.
    - `tpm_unstranded`: Transcripts per million (unstranded).
    - `fpkm_unstranded`: Fragments per kilobase million (unstranded).
    - `fpkm_uq_unstranded`: FPKM upper-quartile normalized (unstranded).

## Usage Notes
- This is a raw data repository. Use these TSV files for downstream bioinformatics tasks such as differential gene expression analysis, clustering, or survival analysis.
- When merging these datasets, ensure you handle the `gene_id` consistently, as each file contains the same set of genes based on the GENCODE reference.
- Some samples may have associated clinical annotation files (`annotations.txt`) which can be mapped back to the case or sample identifiers in the count files.
