#!/usr/bin/env python2.7
#Time-stamp:<Tarela>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import sys
from optparse import OptionParser
import twobitreader
from bx.bbi.bigwig_file import BigWigFile
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def rev(seq):
    revseq = ""
    r = ''
    for i in seq[::-1]:
        if i == 'A':
            r = 'T'
        elif i == 'T':
            r = 'A'
        elif i == 'C':
            r = 'G'
        elif i == 'G':
            r = 'C'
        else:
            print(i)
        revseq += r
    return revseq
def count_cut_nmers( fp, w_plus,w_minus,lflank, rflank ,single_nmer_cutoff,sequence):
    """
    count the number of cuts associated with each nmer in sequence covered by X.
    offset is the position of the cut to be associated with each nmer.
    if offset = 0 the first base of the tag is lined up with the nmer start
    """
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))

    genome = twobitreader.TwoBitFile(sequence)
    # keep count of the number of occurrences of each n-mer
    
    seq_nmer_dict   = {}

    cut_nmer_dict    = {}

    
    for line in fp.readlines():
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        seq = genome[chrm][(start-lflank):(end+rflank)].upper()
        cp = list(w_plus_H.summarize(ll[0],start,end,end-start).sum_data)
        cn = list(w_minus_H.summarize(ll[0],start,end,end-start).sum_data)
        #each = (len(ll)-5)/2
        #cp = (map(float,ll[5:(5+each)]))
        #cn = (map(float,ll[(5+each):(5+each*2)]))

        for k in range( len(cp) ):
            
            p_cut = cp[k]
            n_cut = cn[k]
            
            p_seq = seq[k:(k+lflank+rflank)]
            n_seq = seq[(k+1):(k+lflank+rflank+1)]
       #     rev_n_seq = rev(n_seq)
            if 'N' not in p_seq and p_cut <= single_nmer_cutoff :
                try: 
                    cut_nmer_dict[ p_seq ] += p_cut
                except:
                    cut_nmer_dict[ p_seq ]  = p_cut
                try: 
                    seq_nmer_dict[ p_seq ] += 1
                except:
                    seq_nmer_dict[ p_seq ]  = 1
            if 'N' not in n_seq and n_cut <= single_nmer_cutoff:
                rev_n_seq = rev(n_seq)
                try:
                    cut_nmer_dict[ rev_n_seq ] += n_cut
                except:
                    cut_nmer_dict[ rev_n_seq ]  = n_cut
                try: 
                    seq_nmer_dict[ rev_n_seq ] += 1
                except:
                    seq_nmer_dict[ rev_n_seq ]  = 1
    return seq_nmer_dict, cut_nmer_dict
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
        nmer_seq[seq] = -1
    del allseq
    return nmer_seq


def build_nmer_model(fp, w_plus,w_minus,out,single_nmer_cut,sequence,lflank,rflank):
    """
    count cuts within each nmer sequence 
    """
    # |0.1.2.3.4 lflank=0 rflank=5
    #  0|1.2.3.4 lflank=1 rflank=4
    nmer = lflank + rflank

    # count  
    seq_nmer_dict, cut_nmer_dict = count_cut_nmers( fp, w_plus,w_minus, lflank,rflank ,single_nmer_cut,sequence)
    # -1 : no seq , 0 : seq no cut
    dall = make_nmer_dict(nmer)
    #print len(allcut_nmer_dict),len(cut_nmer_dict_p),len(cut_nmer_dict_n)
    # normalize
    for elem in seq_nmer_dict:
    
        if not dall.has_key(elem):
            continue
        try : 
            val = cut_nmer_dict[elem]
        except:
            val = 0

        
        dall[elem]  =  val / (1.0*seq_nmer_dict[elem]) 
    outf = open(out,'w')
    for seqtype in sorted(dall.keys()):
        pout = str(dall[seqtype])
        nout = str(dall[rev(seqtype)])
        #pcut = str(cut_nmer_dict_p[seqtype])
        #ncut = str(cut_nmer_dict_n[seqtype])
        #bgseq = str(allcut_nmer_dict[seqtype])
        outf.write("\t".join([seqtype,pout,nout])+"\n")
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
    optparser.add_option("-i","--input_peak",dest="peak",type="str",
                         help="peak/interval file in bed format, the region you want to take to calculate DNase sequence bias")
    optparser.add_option("-o","--output",dest="out",type="str",
                         help="")
    optparser.add_option("-p","--bwplus",dest="wp",type="str",
                         help="cleavage bigWiggle file for DNase positive cut")
    optparser.add_option("-n","--bwminus",dest="wm",type="str",
                         help="cleavage bigWiggle file for DNase negative cut")
    
    optparser.add_option("--genome",dest="genome",type="str",default = '/home/sh430/Data/Genome/hg19.2bit',
                         help="twobit sequence file , default : /home/sh430/Data/Genome/hg19.2bit")
    optparser.add_option("--sc",dest="single_nmer_cut",type="int",default = 20,
                         help="cutoff for cleavage count at single loci , if there's more than N cut at a single loci , we think its caused by amplification ,and skip it")
    optparser.add_option("--left",dest="left",type="int",default = 3,
                         help="region length(left) for calculate seqbias , default =3  nmer = left+right")
    optparser.add_option("--right",dest="right",type="int",default = 3,
                         help="region length(right) for calculate seqbias , default =3  nmer = left+right")

                         
#========minor options=============

    (options,args) = optparser.parse_args()

    peak = options.peak
    out = options.out
    w_plus = options.wp
    w_minus = options.wm
    genome = options.genome
    single_nmer_cut = options.single_nmer_cut

    if not peak:
        optparser.print_help()
        sys.exit(1)
    
    fp = open(peak)
    build_nmer_model(fp, w_plus,w_minus,out,single_nmer_cut,genome,options.left,options.right)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



