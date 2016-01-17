# NPG-explorer, Nucleotide PanGenome explorer
# Copyright (C) 2012-2016 Boris Nagaev
#
# See the LICENSE file for terms of use.

# This script reads standard input, writes to standard output

import sys

previous_fields = None

for line in sys.stdin:
    line = line.strip()
    fields = line.split('\t')
    # recover values of "dotted" fields from previous line
    for i in range(4):
        if fields[i] == '.':
            fields[i] = previous_fields[i]
    block = fields[0]
    fragment = fields[1]
    start = fields[2]
    stop_or_letter = fields[3]
    if fragment == 'fragment':
        # skip header line
        continue
    start = int(start)
    format_line = "%(block)s %(fragment)s %(pos)d %(change)s"
    try:
        stop = int(stop_or_letter)
        # gap
        for pos in range(start, stop + 1):
            print(format_line % {
                'block': block,
                'fragment': fragment,
                'pos': pos,
                'change': '-'})
    except:
        # substitution or gap of length 1
        new_letter = stop_or_letter
        assert len(new_letter) == 1
        print(format_line % {
            'block': block,
            'fragment': fragment,
            'pos': start,
            'change': new_letter})
    previous_fields = fields
