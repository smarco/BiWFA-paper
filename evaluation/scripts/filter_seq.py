# Usage:
#    cat *.seq | python3 scripts/filter_seq.py 0 10000

import sys

min_len = int(sys.argv[1])
max_len = int(sys.argv[2])

for line in sys.stdin:
    line1 = line.strip()
    line2 = sys.stdin.readline().strip()

    seq1 = line1[1:]
    seq2 = line2[1:]
    if min_len <= len(seq1) <= max_len and min_len <= len(seq2) <= max_len:
        print(line1)
        print(line2)
