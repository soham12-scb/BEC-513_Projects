#!/usr/bin/env python3
import sys
import re

# Axis boundaries
PM, PX = -500, 499
LM, LX = 50, 400

positions = range(PM, PX + 1)
lengths = range(LM, LX + 1)

pos_count = len(positions)
matrix = [[0] * pos_count for _ in lengths]

# Regex pattern to extract genomic coordinates from column 5 (e.g., "chr1:1000-2000")
region_pattern = re.compile(r"^([^:|]+):(\d+)-(\d+)")

# Caching identical consecutive motif sites to eliminate redundant regex matches
last_site = ""
site_center = 0

for line in sys.stdin:
    columns = line.split()

    # Need at least 10 columns to access column 9/10 (indices 8 and 9)
    if len(columns) < 10:
        continue

    # 1. Landmark motif center logic derived from vplot_col5
    raw_site = columns[4]
    
    if raw_site != last_site:
        match = region_pattern.match(raw_site)
        if not match:
            continue
            
        last_site = raw_site
        ref_start = int(match.group(2))
        ref_end = int(match.group(3))
        
        # Calculate reference midpoint
        site_center = (ref_start + ref_end) // 2

    # 2. Fragment start and end from columns 9 and 10 (indices 8 and 9)
    try:
        start = int(columns[8])
        end = int(columns[9])
    except ValueError:
        continue

    # Calculate properties
    length = end - start
    position = (start + end) // 2 - site_center

    # 3. Populate matrix if within boundaries
    if PM <= position <= PX and LM <= length <= LX:
        matrix[length - LM][position - PM] += 1

# Buffered stream output (faster than repeated print calls)
write = sys.stdout.write
write("length\t" + "\t".join(map(str, positions)) + "\n")

for i, length in enumerate(lengths):
    row_data = "\t".join(str(c) for c in matrix[i])
    write(f"{length}\t{row_data}\n")