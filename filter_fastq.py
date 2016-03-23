#!/usr/bin/env python

# filter_fastq.py

# filter sequence ids by identifiers in a file.

# SNF March 2016

import sys, argparse, re, gzip, os, threading, time
from itertools import islice, izip
import multiprocessing

parser = argparse.ArgumentParser(description="Filter fastq file") 

parser.add_argument("-i", "--input", help="File to filter unpaired reads from (default stdin)", default="stdin") 
parser.add_argument("-1", "--read1", help="Paired-end fastq file one")
parser.add_argument("-2", "--read2", help="Paired-end fastq file two") 
parser.add_argument("-p", "--num-threads", help="Number of threads to start (currently for PE reads only)", type=int)
parser.add_argument("-o", "--output", help="Output filename (default stdout; interpreted as a basename for PE reads)", default="stdout") 
parser.add_argument("-f", "--filter_file", help="File containing strings to include based on IDs in fastq") 
parser.add_argument("-v", "--invert", help="Invert match (exclude matching entries)", action="store_true") 
parser.add_argument("--gzip", help="gzip compress the output", action="store_true") 

args = parser.parse_args() 

if args.input != "stdin" and (args.read1 or args.read2):
    sys.exit("Fatal: please provide only one of --input or --read1/--read2 but not both") 

if (args.read1 and not args.read2) or (args.read2 and not args.read1): 
    sys.exit("Fatal: both --read1 and --read2 must be provided when processing paired-end reads")

if args.read1: # reads are paired now if this is true
    paired = True
else: 
    paired = False

if not paired: 
    if (args.output == "stdout" or args.output == "-"): 
        outfile = sys.stdout
    else: 
        if args.gzip: 
            outfile = gzip.open(args.output, "wb") 
        else:
            outfile = open(args.output, "w") 

    if (args.input == "stdin" or args.input == "-"): 
        infile = sys.stdin
    else: 
        filename, file_extension = os.path.splitext(args.input)
        if file_extension == ".gz":
            infile = gzip.open(args.input, "rb")
        elif file_extension == ".fastq":
            infile = open(args.input, "r") 
        else:
            sys.exit("Fatal: unknown file extension for filename %s" % args.input) 

else: # paired
    filename, file_extension = os.path.splitext(args.read1)
    if file_extension == ".gz":
        infile_read1 = gzip.open(args.read1, "rb") 
    elif file_extension == ".fastq": 
        infile_read1 = open(args.read1, "r") 
    else:
        sys.exit("Fatal: unknown file extension for filename %s" % args.read1) 

    filename, file_extension = os.path.splitext(args.read2)
    if file_extension == ".gz":
        infile_read2 = gzip.open(args.read2, "rb") 
    elif file_extension == ".fastq":
        infile_read2 = open(args.read2, "r") 
    else:
        sys.exit("Fatal: unknown file extension for filename %s" % args.read2) 
    
    if args.gzip: 
        outfile_read1 = gzip.open(args.output + ".1.fastq.gz", "wb") 
        outfile_read2 = gzip.open(args.output + ".2.fastq.gz", "wb") 
    else:
        outfile_read1 = open(args.output + ".1.fastq", "w") 
        outfile_read2 = open(args.output + ".2.fastq", "w") 

    sys.stderr.write("Filtering reads from paired files %s and %s.\n" % (args.read1, args.read2)) 
    
filter_list = [i.strip() for i in open(args.filter_file)]

sys.stderr.write("Read %d identifiers to filter.\n" % len(filter_list))

# note this function uses filter_list as defined globally earlier. can't pass it in b/c need one argument for pool.imap
def filter_pair(lines):
    # lines is a set of four lines from both files that are stitched together. Order is: 
    #   file1line1, file2line1, file1line2, file2line2, etc...
    seqid1 = lines[0].split()[0]
    seqid2 = lines[1].split()[0]

    if seqid1 != seqid2:
        print "WARNING: out of sync pairs detected with ids %s %s" % (seqid1, seqid2)

    if args.invert: 
        if seqid1 in filter_list: 
            #print "Filtering ID r1 %s" % seqid1
            return False
        else:
            return lines
    else:
        if seqid1 in filter_list: 
            #print "Keeping ID r1 %s" % seqid1
            return lines
        else: 
            return False

def filter_read(lines): 
    # lines is a set of four lines from the file.
    seqid = lines[0].split()[0]

    if args.invert: 
        if seqid in filter_list: 
            #print "Filtering ID r1 %s" % seqid1
            return False
        else:
            return lines
    else:
        if seqid in filter_list: 
            #print "Keeping ID r1 %s" % seqid1
            return lines
        else: 
            return False

n_reads = 0 
n_filtered = 0 
pool = multiprocessing.Pool(processes=args.num_threads) 

if paired: 
    # the izip below makes an iterator that loops over both files in sets of four. should be parallel-safe
    res = pool.imap(filter_pair, izip(*[infile_read1, infile_read2]*4))
else: 
    # one-file iterator 
    res = pool.imap(filter_read, izip(*[infile]*4))

while 1: 
    if not n_reads % 1000000:
        sys.stderr.write("Processed %d reads.\n" % n_reads)
    try: 
        lines = res.next() 
        n_reads += 1
        if lines: 
            if paired: 
                outfile_read1.write(lines[0] + lines[2] + lines[4] + lines[6])
                outfile_read2.write(lines[1] + lines[3] + lines[5] + lines[7])
            else: 
                outfile.write("".join(lines))
            if not args.invert: 
                n_filtered += 1
        else: # don't output this
            if args.invert: 
                n_filtered += 1


    except StopIteration: 
        break

pool.close() 
pool.join() 

sys.stderr.write("Processed %d reads, filtered %d reads. (%.2f%%).\n" % (n_reads, n_filtered, n_filtered/float(n_reads)))

