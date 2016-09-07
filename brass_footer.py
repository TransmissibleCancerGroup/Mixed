#!/usr/bin/env python

DOC='''
This program extracts the footer information from BRASS output files.
The output is written as tab-separated fields to the system standard output.

Usage:
    brass_footer.py FILE [FILES...]  # output written to screen
    brass_footer.py FILE [FILES...] > OUTPUT.TSV  # output written to OUTPUT.TSV
    ls PATH/*.BRASSOUT | brass_footer.py # input through pipe

'''
import sys
import os
import logging
from functools import reduce

logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())


class Result(object):
    def __init__(self, filename):
        self.values = {'FileName': filename}

    def add(self, item, value):
        self.values[item] = value

    def keys(self):
        return set(self.values)

    def output(self, keys):
        fields = []
        for k in keys:
            fields.append(str(self.values.get(k, 0)))
        return '\t'.join(fields) + '\n'


class SetOfResults(object):
    def __init__(self):
        self.results = []

    def __len__(self):
        return len(self.results)

    def __str__(self):
        return self.output()

    def keys(self):
        return reduce(set.union, [result.keys() for result in self.results], set())

    def output(self):
        lines = []
        keys = ['FileName'] + [key for key in sorted(self.keys()) if key != 'FileName']
        if keys == ['FileName']:
            raise ValueError('No files had any footer information')
        lines.append('\t'.join(keys) + '\n')
        for result in self.results:
            lines.append(result.output(keys))
        return ''.join(lines)

    def append(self, result):
        self.results.append(result)


def skip_to_end(f):
    """
    Skip towards the end of the file.

    Brass output files can be big, but the footer is always around the
    same size. This makes it faster to skip to the end of the file,
    then rewind to an area just before the footer (about 2000 bytes
    seems to work).
    """
    f.seek(0, os.SEEK_END)
    here = f.tell()
    f.seek(max(0, here-2000), os.SEEK_SET)


def process_file(filename):
    import re

    # Regex extracts FIELD NAME and VALUE from lines of the form
    # "#   string FIELD NAME: int VALUE"
    rgx = re.compile('^[#|%]\s+(.+)\:\s+(\d+)') # comment char can be # or %
    if not os.path.exists(filename):
        logger.error("File {} does not exist".format(filename))
        return

    result = Result(filename)

    with open(filename) as fl:
        skip_to_end(fl)
        for line in fl:
            search = rgx.search(line)
            if search is not None:

                # Remove spaces and put in title case
                itemname = re.sub('[^a-zA-Z<]+', '', search.group(1).title())

                # One field starts with non-alphanumeric chars - " < N read pairs"
                # so this replaces "< N" with TooFew
                itemname = re.sub('<\d?', 'TooFew', itemname)

                # Value is an integer
                value = int(search.group(2))
                result.add(itemname, value)
    return result


if __name__ == "__main__":
    processed = SetOfResults()

    if len(sys.argv) > 1:  # Prefer files on command line
        filenames = sys.argv[1:]

    else:  # No files on the command line, fallback check stdin
        if sys.stdin.isatty():  # Waiting for terminal input, so write some help
            sys.stderr.write(DOC)
            sys.stderr.write("Enter files to scan (CTRL-D when done, CTRL-C to quit):\n")

        try:  # Grab filenames from stdin
            filenames = [f.strip() for f in sys.stdin.read().split() if f.strip() > '']

        except KeyboardInterrupt: # Exit cleanly from CTRL-C, no backtrace
            sys.exit()

    for filename in filenames:
        result = process_file(filename)
        if result is not None:
            processed.append(result)

    if len(processed) > 0:  # Successfully processed some files, so write to stdout
        try:
            sys.stdout.write(str(processed))
        except ValueError as err:
            logger.error(str(err))
            sys.stderr.write(DOC)

    else:  # FAIL!
        logger.error("No files were found")
        sys.stderr.write(DOC)

