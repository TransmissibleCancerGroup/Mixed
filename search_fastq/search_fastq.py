#!/usr/bin/env python
import os, re, sys

def readfq(fp): # this is a generator function
    """
    Function comes from Heng Li: 
    https://raw.githubusercontent.com/lh3/readfq/master/readfq.py

    parameter fp = open file object
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, help='Fastq file to search')
    parser.add_argument('regex', type=str, help='Pattern to search for in the record\'s sequence')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    
    if not os.path.exists(args.filename):
        sys.stderr.write('File {} not found: exiting\n'.format(args.filename))
        sys.exit()

    try:
        regex = re.compile(args.regex)
    except:
        sys.stderr("Error compiling this pattern: '{}' into a regular expression\n".format(args.regex))
        sys.exit()

    with open(args.filename) as fl:
        try:
            for name, seq, qual in readfq(fl):
                if regex.search(seq):

                    # This prints the full fastq record for each match
                    print('@{}\n{}\n+\n{}'.format(name, seq, qual))
        except:
            sys.stderr.write('There was an error processing {} as a FastQ file\n'.format(args.filename))
            sys.stderr.write('The error occurred at file position {} (bytes)\n'.format(fl.tell))
            sys.exit()


