# Project 1: BED File Chromosome Prioritization and Sorting Tool

## Overview
**Project 1** is a high-performance Python utility designed for processing genomic data in standard BED format. It filters and re-sorts BED records according to a custom chromosome priority ranking list (provided via `selection.tsv`). 

The script is optimized for streaming large-scale genomic datasets efficiently via standard input (`stdin`) and standard output (`stdout`), making it suitable for bioinformatics pipelines.

---

## Project Structure
```text
Project 1/
│
├── Final Code P-1.py                       # Main processing script
├── selection.tsv                           # Priority ranking list of chromosomes/regions (input required)
└── sorted_bed_file_per_sample/
    └── sampleX_standard.sorted.bed         # Example sorted BED input file for testing
```

---

## Key Features & Optimizations
1. **Streaming Architecture**: Reads input line-by-line from `stdin` and writes sorted output directly to `stdout`, minimizing memory overhead and enabling UNIX pipe integration.
2. **Fast Lookups**: Loads priority ranks into an in-memory dictionary (`chrom_order`) and caches methods (`dict.get`, `list.append`) to maximize execution speed during inner-loop parsing.
3. **Optimized Parsing**: Splits only the first 3 columns of each BED record, leaving columns 4 onwards untouched to preserve complex metadata, scores, or names.
4. **Custom Sorting**: Sorts filtered records in-place based on user-defined chromosomal priorities specified in `selection.tsv`.

---

## Prerequisites
- **Python 3.x** (No external third-party libraries required; uses standard Python libraries).

---

## Input Formats

### 1. `selection.tsv`
A plain-text file listing target chromosomes or regions in the exact order of desired priority (one identifier per line):
```text
chr1
chr2
chrX
...
```

### 2. BED Input (via `stdin`)
Standard 3-to-6+ column BED format (tab-separated):
```text
chr1    10000   10061   61      .       +
chr1    10000   10086   A00738:475:HCY7JDSX5:4:2669:7256:18568_AGGCCTGTTAAT    86      +
```

---

## Usage

You can run the script from the command line using standard UNIX/PowerShell redirection and pipes:

### In PowerShell (Windows):
```powershell
Get-Content "sorted_bed_file_per_sample\sampleX_standard.sorted.bed" | python "Final Code P-1.py" > "output.sorted.bed"
```

### In Bash / Linux / macOS:
```bash
cat sorted_bed_file_per_sample/sampleX_standard.sorted.bed | python3 "Final Code P-1.py" > output.sorted.bed
```

---

## How It Works
1. **Load Priorities**: Reads `selection.tsv` and assigns a zero-based integer rank to each chromosome.
2. **Stream & Filter**: Iterates over each line from `stdin`, splitting columns to extract the chromosome name. If the chromosome exists in `selection.tsv`, its rank and line are stored as a tuple `(chrom_rank, line)`.
3. **Sort**: Sorts the collected records in-place using Python's efficient Timsort algorithm (primary key: `chrom_rank`).
4. **Output**: Writes the sorted and filtered records sequentially to `stdout`.