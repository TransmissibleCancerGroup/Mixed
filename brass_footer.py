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

logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())

class Result(object):

    def __init__(self, filename):
        self.filename = filename
        self.total_scanned = 0
        self.properly_paired = 0
        self.half_unmapped = 0
        self.near_mate = 0
        self.low_quality = 0
        self.small_insertion = 0
        self.repeat_features = 0
        self.total_groups_found = 0
        self.rearr_groups_omitted = 0
        self.total_groups_emitted = 0

    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            self.filename, self.total_scanned, self.properly_paired,
            self.half_unmapped, self.near_mate, self.low_quality,
            self.small_insertion, self.repeat_features,
            self.total_groups_found, self.rearr_groups_omitted,
            self.total_groups_emitted)


HEADER = '\t'.join([
    "FILENAME",
    "TOTAL_SCANNED",
    "PROPERLY_PAIRED",
    "HALF_UNMAPPED",
    "NEAR_MATE",
    "LOW_QUALITY",
    "SMALL_INSERTION",
    "REPEAT_FEATURES",
    "TOTAL_GROUPS_FOUND",
    "REARR_GROUPS_OMITTED",
    "TOTAL_GROUPS_EMITTED\n"
    ])

def get_num_from_line(line):
    """
    Read the number in the last field of a line, or fail with some
    error info. Failure to read returns 0.
    """
    try:
        num = line.split()[-1]
    except IndexError as err:
        logger.error("Empty line? Line:{} Err:{}".format(line, err))
        return 0

    try:
        num = int(num)
    except ValueError as err:
        logger.error("Not a number? Line:{} Err:{}".format(line, err))
        return 0

    return num

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
    """
    Read pertinent info from the file footer.
    """
    if not os.path.exists(filename):
        logger.error("File {} does not exist".format(filename))
        return

    result = Result(filename)

    with open(filename) as fl:
        skip_to_end(fl)
        for line in fl:
            trim = line.lstrip('#').strip()
            if trim == '':
                continue

            if trim.startswith('Total reads scanned'):
                result.total_scanned = get_num_from_line(trim)

            elif trim.startswith("Properly paired"):
                result.properly_paired = get_num_from_line(trim)

            elif trim.startswith("(Half-)unmapped"):
                result.half_unmapped = get_num_from_line(trim)

            elif trim.startswith("Near mate"):
                result.near_mate = get_num_from_line(trim)

            elif trim.startswith("Low quality"):
                result.low_quality = get_num_from_line(trim)

            elif trim.startswith('Small insertion'):
                result.small_insertion = get_num_from_line(trim)

            elif trim.startswith('Repeat features'):
                result.repeat_features = get_num_from_line(trim)

            elif trim.startswith('Total groups found'):
                result.total_groups_found = get_num_from_line(trim)

            elif trim.startswith('< 2 read pairs'):
                result.rearr_groups_omitted = get_num_from_line(trim)

            elif trim.startswith('Total groups emitted'):
                result.total_groups_emitted = get_num_from_line(trim)

    return result


if __name__ == "__main__":
    processed = []

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
        sys.stdout.write(HEADER)
        for result in processed:
            sys.stdout.write(str(result))

    else:  # FAIL!
        logger.error("No files were found")
        sys.stderr.write(DOC)

