# Project 2: Fragment Length Distribution & Rescaling

## Overview
**Project 2** is a bioinformatics toolkit designed for analyzing genomic fragment length distributions and performing distribution-based rescaling (downsampling). It is particularly useful for ATAC-seq, ChIP-seq, or cfDNA analysis where fragment length profiles need to be normalized to match a gold-standard reference.

The tool uses **rejection sampling** to filter fragments from a query dataset so that the resulting distribution matches a specified target distribution.

---

## Project Structure
```text
Project 2/
├── Project 2 code.py                   # Script to calculate fragment length frequencies vs a reference
├── rescale_bed.py                      # Script to downsample BED records via rejection sampling
├── reference.hist                      # Input: Target length distribution (gold-standard)
├── query.bed                           # Input: Raw genomic query data (BED format)
├── plot_three_way.gp                   # Gnuplot script for visualization
└── Output Files:
    ├── fragment_length_comparison.txt  # Comparison table: Query vs Reference frequencies
    ├── rescaled_query.bed              # Output: Downsampled/Normalized BED file
    ├── three_way_comparison.txt        # Combined frequency table for all 3 states
    └── three_way_comparison_distribution.png # Final visualization plot
```

---

## Prerequisites
- **Python 3.x** (Standard library only: `sys`, `os`, `random`)
- **Gnuplot** (Optional: Required for generating the comparison plot)

---

## Workflow Guide

Follow these steps to normalize your fragment length distribution and verify the results:

### 1. Initial Comparison
Calculate the distribution of your raw `query.bed` and compare it to `reference.hist`.
```bash
python "Project 2 code.py" < query.bed > fragment_length_comparison.txt
```
This produces `fragment_length_comparison.txt`, which is required for the next step.

### 2. Rescale (Normalize) the BED File
Use the comparison results to downsample the `query.bed` file. The script will filter reads based on the ratio of reference frequencies to query frequencies.
```bash
python rescale_bed.py < query.bed > rescaled_query.bed
```

### 3. Verify the Rescaled Distribution
Calculate the distribution of the newly created `rescaled_query.bed` to confirm it now matches the reference.
```bash
python "Project 2 code.py" < rescaled_query.bed > rescaled_length_comparison.txt
```

### 4. Consolidate Data for Plotting
Merge the original, rescaled, and reference distributions into a single table for visualization.
```python
# Run this snippet in Python
import sys

def load_freqs(filename, col_idx):
    freqs = {}
    with open(filename) as f:
        f.readline() # skip header
        for line in f:
            p = line.strip().split()
            if len(p) >= col_idx + 1:
                freqs[p[0]] = p[col_idx]
    return freqs

orig = load_freqs("fragment_length_comparison.txt", 1)
ref = load_freqs("fragment_length_comparison.txt", 2)
resc = load_freqs("rescaled_length_comparison.txt", 1)

all_lens = sorted(set(orig.keys()) | set(ref.keys()) | set(resc.keys()), key=int)

with open("three_way_comparison.txt", "w") as out:
    out.write("length\torig_freq\tresc_freq\tref_freq\n")
    for l in all_lens:
        out.write(f"{l}\t{orig.get(l, '0')}\t{resc.get(l, '0')}\t{ref.get(l, '0')}\n")
```

### 5. Generate Comparison Plot
Use Gnuplot to create the final 3-way comparison graph.
```bash
gnuplot plot_three_way.gp
```
This will generate `three_way_comparison_distribution.png`.

---

## File Formats

### `reference.hist`
A two-column, tab-separated file containing fragment length and its corresponding probability (normalized to 1).
```text
#Length  Probability
50      0.0012
51      0.0015
...
```

### `query.bed`
Standard BED format. The script uses columns 2 (Start) and 3 (End) to calculate fragment length as `End - Start`.
```text
chr1    10000   10150
chr1    20000   20145
```
