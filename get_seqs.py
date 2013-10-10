#!/usr/bin/python

import sys
import re
import argparse
import urllib2

URL = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=%(db)s&id=%(id)s'+\
      '&format=%(format)s&style=raw'

def get_seqs(args):
    name2name = {}
    format = 'fasta' if (args.type == 'fasta') else 'default'
    for line in args.table:
        line = line.strip()
        if line:
            try:
                fasta_id, genome, chromosone, circular = line.split()
                db = 'refseqn' if re.match(r'^\w\w_', fasta_id) else 'embl'
                url = URL % {'db': db, 'id': fasta_id, 'format': format}
                print(url)
                args.out.write(urllib2.urlopen(url).read())
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
    p.add_argument('--type', help='Type of content downloaded (fasta|genes)',
            default='fasta', choices=('fasta', 'genes'))
    args = p.parse_args()
    get_seqs(args)

if __name__ == '__main__':
    try:
        main()
    except Exception, e:
        print e
        sys.exit(255)

