#!/usr/bin/python

import sys
import re
import argparse
import urllib2

URL = 'http://www.ebi.ac.uk/ena/data/view/%s&display=fasta&expanded=true'

def get_seqs(args):
    name2name = {}
    for line in args.table:
        line = line.strip()
        if line:
            try:
                fasta_id, genome, chromosone, circular = line.split()
                args.out.write(urllib2.urlopen(URL % fasta_id).read())
            except:
                print 'Warning: bad line in table: ' + line

def main():
    r = argparse.FileType('r')
    w = argparse.FileType('w')
    p = argparse.ArgumentParser(
        description='Download sequences and write to fasta file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--table',
            help='Table (fasta_id genome chromosone c[ircular]/l[inear])',
            metavar='FILE', type=r, required=True)
    p.add_argument('--out', help='Output new fasta file',
            metavar='FILE', type=w, required=True)
    args = p.parse_args()
    get_seqs(args)

if __name__ == '__main__':
    try:
        main()
    except Exception, e:
        print e
        sys.exit(255)

