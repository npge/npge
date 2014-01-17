#!/usr/bin/python

import sys

table = sys.argv[1]
block = sys.argv[2]

dist = {}
frs = set()

for line in open(table):
    bl, f1, f2, dst = line.strip().split()
    if bl == block:
        dist[f1, f2] = float(dst)
        dist[f2, f1] = float(dst)
        dist[f1, f1] = 0.0
        dist[f2, f2] = 0.0
        frs.add(f1)
        frs.add(f2)

print("%5i" % len(frs))
for f1 in frs:
    d = ["%10f" % dist[f1, f2] for f2 in frs]
    print("%10s" % f1[:10] + "".join(d))

