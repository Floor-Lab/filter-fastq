# filter-fastq
Filter reads from a FASTQ file using a list of identifiers.

Each entry in the input FASTQ file (or files) is checked against all entries in the identifier list. Matches are included by default, or excluded if the --invert flag is supplied. Paired-end files are kept consistent (in order).

This is almost certainly not the most efficient way to implement this filtering procedure. I tested a few different strategies and this one seemed the fastest. Current timing with 16 processes is about 10 minutes per 1M paired reads with gzip'd input and output, depending on the length of the identifier list to filter by.

# usage 

```
usage: filter_fastq.py [-h] [-i INPUT] [-1 READ1] [-2 READ2] [-p NUM_THREADS]
                       [-o OUTPUT] [-f FILTER_FILE] [-v] [--gzip]

Filter fastq file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        File to filter unpaired reads from (default stdin)
  -1 READ1, --read1 READ1
                        Paired-end fastq file one
  -2 READ2, --read2 READ2
                        Paired-end fastq file two
  -p NUM_THREADS, --num-threads NUM_THREADS
                        Number of threads to start (currently for PE reads
                        only)
  -o OUTPUT, --output OUTPUT
                        Output filename (default stdout; interpreted as a
                        basename for PE reads)
  -f FILTER_FILE, --filter_file FILTER_FILE
                        File containing strings to include based on IDs in
                        fastq
  -v, --invert          Invert match (exclude matching entries)
  --gzip                gzip compress the output
```
# File formats

Currently this only works on FASTQ files (https://en.wikipedia.org/wiki/FASTQ_format).

gzip'd input is accepted and automatically detected.

The identifier list should be a newline-separated list of each identifier to match against. For example: 

```
@K00188:86:H3NM2BBXX:4:1101:1316:1103
@K00188:86:H3NM2BBXX:4:1101:3934:1103
```
