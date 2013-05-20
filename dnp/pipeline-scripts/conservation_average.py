#!/opt/bin/python
# Time-stamp: <2011-07-15 21:12:20 Jian Ma>

"""
Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
@modified by Qin Qian: to extract Phastcon score only without images
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import re
import logging
import subprocess
import math
from optparse import OptionParser

from CistromeAP.taolib.CoreLib.Parser import WiggleIO, BedIO
from CistromeAP.taolib.CoreLib.BasicStat.Func import * 

try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )
#bigWigSummary = 'bigWigSummary'

# ------------------------------------
# Misc functions
# ------------------------------------

error  = logging.critical		# function alias
warn   = logging.warning
debug  = logging.debug
info   = logging.info

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog <-d path> [options] <bed files> ..."
    description = "Draw conservation plot for many bed files."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option('-w',dest='w',type='int',default=1000, help="window width centered at middle of bed regions,default: 1000")
    optparser.add_option('-d','--phasdb',dest='phasdb',help= 'The directory to store phastcons scores in the server')
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    
    (options,bedfiles) = optparser.parse_args()
    options.pf_res = options.w / 100 # get 100 points to plot
    options.w = options.pf_res * 100 # trim

    bedfiles = map(os.path.abspath,bedfiles)
    bedfilenames = map(os.path.basename,bedfiles)

    bedfilenum = len(bedfiles)

    if bedfilenum < 1 or not options.phasdb:
        optparser.print_help()
        sys.exit(1)

    # check the files
    for f in bedfiles:
        if not os.path.isfile(f):
            error("%s is not valid!" % f)
            sys.exit(1)

    # check phastcons db
    if not os.path.isdir(options.phasdb):
        error("%s is not valid!" % options.phasdb)
        sys.exit(1)

    # change wd to phastcons db path
    olddir = os.path.abspath('.')
    os.chdir(options.phasdb)
    phas_chrnames = []
    files_phasdb = os.listdir('.')
    for file_phasdb in files_phasdb:
        if file_phasdb.endswith('.bw'):
            name = file_phasdb.rstrip('.bw')
            phas_chrnames.append(name)

    if not phas_chrnames:
        error("%s has no valid phastcons db bw files!" % options.phasdb)
        sys.exit(1)
        
    info("number of bed files: %d" % bedfilenum)

    avgValues = []

    # for each bed file
    for f in bedfiles:
        info("extract phastcons scores using %s" % f)
        scores = extract_phastcons(f,phas_chrnames, options.w, options.pf_res)
        avgValues.append(scores)
    print(sum(avgValues[0]) / len(avgValues[0]))

def extract_phastcons ( bedfile, phas_chrnames, width, pf_res ):
    """Extract phastcons scores from a bed file.

    Return the average scores
    """
    info("read bed file...")
    bfhd = open(bedfile)
    bed = BedIO.parse_BED(bfhd)

    # calculate the middle point of bed regions then extend left and right by 1/2 width
    bchrs = bed.peaks.keys()
    bchrs.sort()

    chrs = []
    for c in phas_chrnames:
        if c in bchrs:
            chrs.append(c)

    sumscores = []
    for chrom in chrs:
        info("processing chromosome: %s" %chrom)
        pchrom = bed.peaks[chrom]
        bw = BigWigFile(open(chrom+'.bw', 'rb'))
        for i in range(len(pchrom)):
            mid = int((pchrom[i][0]+pchrom[i][1])/2)
            left = int(mid - width/2)
            right = int(mid + width/2)
            
            if left < 0:
                left = 0
                right = width
            
            summarize = bw.summarize(chrom, left, right, width/pf_res)
            if not summarize:
                continue
            dat = summarize.sum_data / summarize.valid_count
            #dat = dat.strip().split('\t')
            sumscores.append(dat)
            
    sumscores = map(list, zip(*sumscores))
    sumscores = [[t2 for t2 in t if not math.isnan(t2)] for t in sumscores]
    try:
        conscores = [sum(t)/len(t) for t in sumscores]
    except ZeroDivisionError:
        conscores = [0] * (width/pf_res)
    return  conscores
        
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
