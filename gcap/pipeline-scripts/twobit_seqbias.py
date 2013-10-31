#Time-stamp:<Tarela>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging
import string
import copy,time
import twobitreader
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def make_nmer_dict(n):
    nmer_seq = {}
    bp = ['A','C','G','T']
    allseq = [0]*n
    allseq[0] = bp
    i=1
    while i < n:
        allseq[i] = []
        for previous_seq in allseq[i-1]:
            for add_bp in bp:
                new_seq = previous_seq + add_bp
                allseq[i].append(new_seq)
        i += 1
    for seq in allseq[n-1]:
  #      print seq
        nmer_seq[seq] = 0
    del allseq
    return nmer_seq
    
def seqbias(peak,tag,out,sequence,flank):
    genome = twobitreader.TwoBitFile(sequence) 
    pcut = make_nmer_dict(2*flank)
    ncut = make_nmer_dict(2*flank)
    bgseq = make_nmer_dict(2*flank)
    inf = open(peak)
    for line in inf:
        ll = line.strip().split("\t")
        seq = genome[ll[0]][(int(ll[1])-flank):(int(ll[2])+flank)]
        for i in range(len(seq)-2*flank):
            subseq = seq[i:(i+2*flank)]
            if bgseq.has_key(subseq):
                bgseq[subseq]+=1
            else:
                pass
    inf.close()
    inf = open(tag)
    for line in inf:
        ll = line.strip().split("\t")
        if len(ll) <6:
            continue
        if ll[5] == '+':
            seq = genome[ll[0]][(int(ll[1])-3):(int(ll[1])+3)].upper()
            if pcut.has_key(seq):
                pcut[seq]+=1
            else:
                ## print(seq) ## for reads shorter than 6
                pass
        elif ll[5] == '-':
            seq = genome[ll[0]][(int(ll[2])-1-2):(int(ll[2])-1+4)].upper()
            if ncut.has_key(seq):
                ncut[seq]+=1
            else:
                ## print(seq) ## for reads shorter than 6
                pass
    inf.close()
    outf = open(out,'w')
    for seqtype in sorted(pcut.keys()):
        pbias = float(pcut[seqtype])/float(bgseq[seqtype])
        nbias = float(ncut[seqtype])/float(bgseq[seqtype])
        outf.write("\t".join(map(str,[seqtype,pbias,nbias,float(pcut[seqtype]),float(ncut[seqtype]),float(bgseq[seqtype])]))+"\n")
    outf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-p","--peak",dest="peak",type="str",
                         help="")
    optparser.add_option("-t","--tag",dest="tag",type="str",
                         help="")              
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")              
         
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/home/sh430/Data/Genome/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("-f","--flank",dest="flank",type="int",default=3,
                         help="flanking region for n-mer , default =3 means n=6")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    peak = options.peak
    tag = options.tag
    out = options.out
    seq = options.sequence
    flank = options.flank
    if not peak:
        optparser.print_help()
        sys.exit(1)
    
    seqbias(peak,tag,out,seq,flank)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


