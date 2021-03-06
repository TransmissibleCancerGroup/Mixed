#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

/*
Extract summary statistics from samtools depth output

Compile with
gcc -O3 -std=c99 -o depthreader depthreader.c -lz -lm
*/


/* Size of the block of memory to use for reading. */
#define LENGTH 0x1000
#define CUTOFF 10000

void clearBuffer(unsigned char *, int);

/*
Wipe buffer by filling with null bytes
*/
void clearBuffer(unsigned char * buffer, int j) {
    for (int i=0; i <= j; i++) {
        buffer[i] = '\0';
    }
}

int main (int argc, const char** argv)
{
    gzFile * file;
    long int nLines=0;
    int nTabs=0;
    long int sumCoverage=0;
    long int sumSqCoverage=0;
    int maxCoverage=0;
    int j=0;
    unsigned char field[LENGTH];
    long int histogram[CUTOFF+1];  // use this to calculate median
    for (size_t i=0; i < (CUTOFF+1); ++i) {
        histogram[i] = 0;
    }

    if (argc != 2) {
        fprintf(stderr, "Usage: %s file\n", argv[0]);
        return 1;
    }

    file = gzopen (argv[1], "r");
    if (! file) {
        fprintf (stderr, "gzopen of '%s' failed: %s.\n", argv[1],
                 strerror (errno));
            exit (EXIT_FAILURE);
    }

    clearBuffer((unsigned char *)field, LENGTH-1);

    while (1) {
        int err;
        int bytes_read;
        int fieldValue;
        unsigned char buffer[LENGTH];
        bytes_read = gzread (file, buffer, LENGTH - 1);
        buffer[bytes_read] = '\0';

        // examine what's in the buffer
        for (int i=0; i < bytes_read; ++i) {
            if (nTabs == 2) {  // Reached second tab, start recording field contents
                field[j] = buffer[i];
                j++;
            }
            if (buffer[i]=='\n') {
                nLines += 1;
                fieldValue = atoi((char *) field);

                // Keep track of max read coverage value seen
                if (fieldValue > maxCoverage) {
                    maxCoverage = fieldValue;
                }

                // Place value in histogram for median
                if (fieldValue >= CUTOFF) {
                    histogram[CUTOFF] += 1;  // collapse all huge values into last bin
                }
                else {
                    histogram[fieldValue] += 1;
                }
                sumCoverage += fieldValue;
                sumSqCoverage += fieldValue*fieldValue;

                //reset field buffer and tab counter ready for next line
                nTabs = 0;
                clearBuffer((unsigned char *) field, j);
                j = 0;
            }
            else if (buffer[i]=='\t') {
                nTabs += 1;
            }
        }

        if (bytes_read < LENGTH - 1) {
            if (gzeof (file)) {
                break;
            }
            else {
                const char * error_string;
                error_string = gzerror (file, & err);
                if (err) {
                    fprintf (stderr, "Error: %s.\n", error_string);
                    exit (EXIT_FAILURE);
                }
            }
        }
    }
    gzclose (file);

    // Now find median

    double m50 = nLines / 2.0;
    long int runningsum = 0;
    int median = 0;

    for (size_t i=0; i < (CUTOFF+1); ++i) {
        runningsum += histogram[i];
        if (runningsum > m50) {
            median = (int) i;  // should be close enough, will be off by 0.5 in the unlikely event the previous bin has exactly m50 reads running total
            break;
        }
    }

    printf ("Read total of %ld lines\n", nLines);
    printf ("Sum coverage was %ld\n", sumCoverage);
    printf ("SumSq coverage was %ld\n", sumSqCoverage);
    printf ("Max coverage was %d\n", maxCoverage);
    printf ("  mean   = %f\n", (double)sumCoverage/(double)nLines);
    printf ("  s.d.   = %f\n", sqrt(((double)sumSqCoverage - ((double)sumCoverage*(double)sumCoverage/(double)nLines)) / ((double)nLines - 1)));
    printf ("  median = %d\n", median);
    return 0;
}
