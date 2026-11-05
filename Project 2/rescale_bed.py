#!/usr/bin/env python3
import sys
import os
import random

# Fixed seed for reproducibility
random.seed(42)

# Auto-detect existing frequency comparison file
if os.path.exists("fragment_length_comparison.txt"):
    FREQ_FILE = "fragment_length_comparison.txt"
elif os.path.exists("post_sampling_comparison.txt"):
    FREQ_FILE = "post_sampling_comparison.txt"
else:
    sys.stderr.write("Error: Neither 'fragment_length_comparison.txt' nor 'post_sampling_comparison.txt' was found.\n")
    sys.exit(1)

# -------------------------------------------------------------
# Step 1: Read query and reference distributions
# -------------------------------------------------------------
p_query = {}
p_ref = {}
min_len = 10**9
max_len = -1

with open(FREQ_FILE, "r") as f:
    header = f.readline()  # Skip header: length  query_freq  ref_freq
    for line in f:
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        try:
            flen = int(parts[0])
            q_val = float(parts[1])
            r_val = float(parts[2])
        except ValueError:
            continue

        p_query[flen] = q_val
        p_ref[flen] = r_val

        if flen < min_len:
            min_len = flen
        if flen > max_len:
            max_len = flen

# -------------------------------------------------------------
# Step 2: Compute rejection sampling probabilities p_keep(L)
# -------------------------------------------------------------
ratios = {}
max_ratio = 0.0

for flen in range(min_len, max_len + 1):
    q = p_query.get(flen, 0.0)
    r = p_ref.get(flen, 0.0)
    if q > 0.0:
        ratio = r / q
        ratios[flen] = ratio
        if ratio > max_ratio:
            max_ratio = ratio
    else:
        ratios[flen] = 0.0

p_keep = {}
if max_ratio > 0.0:
    for flen in range(min_len, max_len + 1):
        p_keep[flen] = ratios[flen] / max_ratio
else:
    sys.stderr.write("Error: max_ratio is 0. Check your frequency distribution file.\n")
    sys.exit(1)

# -------------------------------------------------------------
# Step 3: Stream BED records from stdin and filter to stdout
# -------------------------------------------------------------
write = sys.stdout.write
input_total = 0
accepted_total = 0

for line in sys.stdin:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        continue

    input_total += 1
    cols = stripped.split("\t")
    if len(cols) < 3:
        cols = stripped.split()
    if len(cols) < 3:
        continue

    try:
        # True fragment length is ALWAYS End - Start
        start = int(cols[1])
        end = int(cols[2])
        flen = end - start
    except ValueError:
        continue

    if min_len <= flen <= max_len:
        prob = p_keep.get(flen, 0.0)
        if random.random() < prob:
            write(line)
            accepted_total += 1

sys.stderr.write(f"[rescale_bed] Input fragments: {input_total} | Accepted: {accepted_total} ({accepted_total / max(input_total, 1) * 100:.2f}%)\n")
