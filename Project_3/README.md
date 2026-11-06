# Project 3: V-Plot Generation

## Overview
**Project 3** is a bioinformatics toolkit designed for generating **V-plots** (fragment length versus distance from a genomic landmark). It is commonly used in ATAC-seq, ChIP-seq, and nucleosome positioning analyses to visualize how genomic fragments (such as transcription factor binding sites or nucleosomes) are distributed around reference landmark centers.

The script processes mapped genomic intervals from standard input (BED format), calculates fragment lengths and relative midpoint positions, and aggregates them into a 2D frequency matrix.

---

## Project Structure
```text
Project 3/
├── Project 3 code.py       # Main script to compute the V-plot matrix
├── mapped.bed              # Input: Mapped genomic fragments (BED format)
├── mapped.bed.gz           # Compressed input fragments (alternative)
├── vplot_matrix.txt        # Output: 2D frequency matrix (Length vs. Relative Position)
└── vplot_gnuplot.png       # Output: Visualization plot of the V-plot matrix
```

---

## Prerequisites
- **Python 3.x** (Standard library only: `sys`, `re`)
- **Gnuplot** (Optional: For rendering the V-plot matrix into an image)

---

## Workflow & Usage

### 1. Generate the V-Plot Matrix
Pipe your mapped BED file into the Python script to generate the 2D frequency matrix (`vplot_matrix.txt`).

```bash
python "Project 3 code.py" < mapped.bed > vplot_matrix.txt
```

*Note:* If using compressed input:
```bash
gunzip -c mapped.bed.gz | python "Project 3 code.py" > vplot_matrix.txt
```

### 2. Matrix Parameters
By default, the script analyzes:
- **Relative Position (X-axis):** `-500` bp to `+499` bp relative to the reference landmark center.
- **Fragment Length (Y-axis):** `50` bp to `400` bp.

---

## Input Format Requirements
The input should be a tab-separated BED-like format containing at least 10 columns:
- **Column 5 (Index 4):** Reference landmark region (e.g., `chr1:1000-2000`) used to calculate the landmark center midpoint.
- **Column 9 & 10 (Indices 8 & 9):** Fragment start and end coordinates, used to compute fragment length (`End - Start`) and midpoint position relative to the landmark.
