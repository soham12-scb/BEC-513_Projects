import sys

# 1. Load priority ranks
chrom_order = {}
rank = 0
with open("selection.tsv", "r") as f:
    for line in f:
        chrom = line.strip()
        if chrom:
            chrom_order[chrom] = rank
            rank += 1

# 2. Local variable lookups for maximum loop speed
records = []
append_record = records.append
order_get = chrom_order.get

# 3. Stream and parse with fast splitting
for line in sys.stdin:
    line = line.strip()
    if not line:
        continue

    # Split only the first 3 columns; keep the rest untouched
    parts = line.split("\t", 3)
    chrom_rank = order_get(parts[0])

    if chrom_rank is not None:
        append_record((chrom_rank, line))

# 4. Sort in-place
records.sort()

# 5. Fast block writing to stdout
write_out = sys.stdout.write
for item in records:
    write_out(item[1] + "\n")