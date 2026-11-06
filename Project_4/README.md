# Project 4: TCGA-LUAD NKX2-1 Expression Analysis

## Overview
**Project 4** is a bioinformatics pipeline focused on analyzing **NKX2-1** gene expression in Lung Adenocarcinoma (TCGA-LUAD) datasets. The workflow parses RNA-Seq count files downloaded from the Genomic Data Commons (GDC), performs statistical testing (Welch's two-sample t-test), and generates visualization plots.

## Project Structure
```text
Project 4/
├── project4 (1).R                      # R script for extraction, statistics, and plotting
├── gdc_sample_sheet.tsv                # Metadata manifest mapping sample types to files
├── README.md                           # Project documentation
└── [GDC Case Folders]/                  # Subdirectories containing STAR count TSVs
    └── *.rna_seq.augmented_star_gene_counts.tsv
```

---

## Prerequisites (Linux)
- **R (v4.0 or higher)**: Standard base R installation (`Rscript` binary in PATH).
- **No extra R packages required**: The script relies exclusively on standard base R libraries (`stats`, `graphics`, `grDevices`).

To verify R is installed on your Linux system:
```bash
Rscript --version
```

---

## Running the Pipeline on Linux

### Step 1: Navigate to Project Directory
```bash
cd /path/to/Project\ 4
```

### Step 2: Ensure Archive Extraction (If applicable)
If your raw GDC data is compressed in a `.tar.gz` archive, extract it first:
```bash
tar -xzf tcga_data.tar.gz
```

### Step 3: Execute the R Script
Run the analysis script via `Rscript`:
```bash
Rscript "Project_4_code"
```

---

## Outputs Generated

Upon execution, the script generates two main output files in the root directory:

1. **`NKX2-1_expression_values.tsv`**
   - Contains extracted raw TPM, `log2(TPM + 1)` values, sample types (*Primary Tumor* vs. *Solid Tissue Normal*), and file paths for all samples.

2. **`NKX2-1_expression.png`**
   - High-resolution boxplot comparing `log2(TPM + 1)` expression across sample groups with an annotated p-value from Welch's t-test.

---

## Additional OS Notes

<details>
<summary>Running on Windows (PowerShell)</summary>

```powershell
cd "C:\Users\soham\Downloads\Project 4"
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "project4 (1).R"
```
</details>
