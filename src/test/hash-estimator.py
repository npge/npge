import sys
s = set()
bad = sum = 0
for line in sys.stdin:
    n = int(line)
    if n in s:
        bad += 1
    sum += 1
    s.add(n)
print(float(bad)/sum)

