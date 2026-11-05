#!/usr/bin/env python3
import sys
import os

REF_FILE = "reference.hist"

if not os.path.exists(REF_FILE):
    sys.stderr.write(f"Error: '{REF_FILE}' not found in current directory.\n")
    sys.exit(1)

# -------------------------------------------------------------
# Step 1: Read reference.hist
# -------------------------------------------------------------
ref_prob = {}
min_len = 10**9
max_len = -1

with open(REF_FILE, "r") as f:
    for line in f:
        parts = line.strip().split()
        if not parts or parts[0].startswith("#"):
            continue
        try:
            flen = int(parts[0])
            prob = float(parts[1])
            ref_prob[flen] = prob
            if flen < min_len:
                min_len = flen
            if flen > max_len:
                max_len = flen
        except ValueError:
            continue

if not ref_prob:
    sys.stderr.write("Error: reference.hist is empty or malformed.\n")
    sys.exit(1)

# -------------------------------------------------------------
# Step 2: Stream input BED from stdin
# -------------------------------------------------------------
num_lengths = max_len - min_len + 1
counts = [0] * num_lengths
total_valid = 0
total_read = 0

for line in sys.stdin:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        continue

    total_read += 1
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
        counts[flen - min_len] += 1
        total_valid += 1

# Progress metric to stderr
sys.stderr.write(f"[Project2t] Lines read: {total_read} | Valid fragments counted: {total_valid}\n")

# -------------------------------------------------------------
# Step 3: Stream output comparison table to stdout
# -------------------------------------------------------------
sys.stdout.write("length\tquery_freq\tref_freq\n")

if total_valid > 0:
    for flen in range(min_len, max_len + 1):
        q_freq = counts[flen - min_len] / total_valid
        r_freq = ref_prob.get(flen, 0.0)
        sys.stdout.write(f"{flen}\t{q_freq:.8e}\t{r_freq:.8e}\n")
else:
    sys.stderr.write("[Project2t] Warning: 0 valid reads were found. No data rows written.\n")
