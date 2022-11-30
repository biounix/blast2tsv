#!/usr/bin/env python3

# blast2tsv.py blast.tsv genome.lengths --evalue e-value --extend length --merge --window length --output output.tsv
# blast2tsv.py blast.tsv genome.lengths --evalue 0.001   --extend 50000  --merge --window 1000   --output output.tsv

import argparse


# Argument parser.
parser = argparse.ArgumentParser(description = 'Parse tabular Blast output file.')

parser.add_argument('blast'          , metavar = 'blast.tsv'     , type = str                          , help = 'blast tabular output'                                                 )
parser.add_argument('lengths'        , metavar = 'genome.lengths', type = str                          , help = 'contig lengths'                                                       )
parser.add_argument('-e', '--evalue' , metavar = 'e-value'       , type = float, default = 1           , help = 'e-value cutoff'                                                       )
parser.add_argument('-o', '--output' , metavar = 'output.tsv'    , type = str  , default = 'output.tsv', help = 'output TSV file'                                                      )
parser.add_argument('-x', '--extend' , metavar = 'length'        , type = int  , default = 0           , help = 'number of bases to extend the sequence before start and end positions')
parser.add_argument('-w', '--window' , metavar = 'window size'   , type = int  , default = 0           , help = 'maximum length between hits to be merged.'                            )

parser.add_argument('-m', '--merge'  , action = argparse.BooleanOptionalAction                         , help = 'merge overlapping hits in a single region')

parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s 1.0.0')

args = parser.parse_args()

# Reading user input file.
cutoff = args.evalue;
extend = args.extend;
merge  = args.merge ;
window = args.window;


lens = {}

# https://stackoverflow.com/a/11198762
with open(args.lengths) as lengths:

    for line in lengths:
        line = line.rstrip() # https://stackoverflow.com/a/275025

        l = line.split()
        # l = line.split('\t')

        length = int(l[0])
        contig =     l[1]

        lens[contig] = length


blastdic = {}

with open(args.blast) as blast:

    for line in blast:
        if line.startswith('#'):
            continue

        line = line.rstrip() # https://stackoverflow.com/a/275025

        # l = line.split()
        l = line.split('\t')

        # protein = l[ 0]
        contig  =       l[ 1]
        start   =   int(l[ 8])
        end     =   int(l[ 9])
        evalue  = float(l[10]) # Check column!

        if evalue > cutoff:
            continue

        # strand = 'fwd'

        # Swap start and end positions for reverse alignments.
        if start > end:
            # https://stackoverflow.com/a/14836456
            start, end = end, start

            # strand = 'rev'

        # https://stackoverflow.com/a/44035382/20592501
        if contig in blastdic:
            # https://stackoverflow.com/q/252703/20592501
            blastdic[contig][0].append(start)
            blastdic[contig][1].append(end  )
            # blastdic[contig][2].append(strand)
        else:
            blastdic[contig] = [ [start], [end] ]
            # blastdic[contig] = [ [start], [end], [strand] ]


def myprint(file, contig, start, end):
    s = max(start - extend, 0           );
    e = min(end   + extend, lens[contig]);

    out = (contig, s, e - s);

    # https://stackoverflow.com/a/3590168
    file.write('\t'.join( map(str, out) ) + '\n')


with open(args.output, 'w') as output:

    for contig in blastdic:
        slist = blastdic[contig][0]
        elist = blastdic[contig][1]

        # https://stackoverflow.com/a/9764364/20592501
        # slist, elist = zip( *sorted( zip(slist, elist) ) )

        # https://stackoverflow.com/a/8372442/20592501
        slist, elist = ( list(a) for a in zip( *sorted( zip(slist, elist) ) ) )
        # list( map(list, zip( *sorted( zip(slist, elist) ) ) ) )

        if merge:
            start = slist.pop(0)
            end   = elist.pop(0)

            # https://stackoverflow.com/a/522578/20592501
            for i, new_start in enumerate(slist):
                new_end = elist[i]

                if new_start > end + window:
                    myprint(output, contig, start, end)

                    start = new_start

                if new_end > end:
                    end = new_end

            myprint(output, contig, start, end)


        else:
            for i, start in enumerate(slist):
                end = elist[i]

                myprint(output, contig, start, end)

