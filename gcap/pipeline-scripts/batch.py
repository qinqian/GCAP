#!/usr/bin/env python

"""
run gcap pararell
"""
from optparse import OptionParser
from ConfigParser import ConfigParser
import os
import sys

def main():
    usage = "usage: %prog -c <conf> -f <db>"
    description = "Draw conservation plot for many bed files."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option('-c',dest='c', help="config template")
    optparser.add_option('-f',dest='f', help="db: id    species treat   control")

    (options, args) = optparser.parse_args()

    if not options.c or not options.f:
        optparser.print_help()
        sys.exit(1)
        
    cf = ConfigParser()
    cf.optionxform = str
    cf.read(options.c)

    with open(options.f) as content:
        for line in content:
            line = line.strip().split()
            cf.set("Basis", "id", line[0])
            cf.set("Basis", "species", line[1])
            cf.set("Basis", "treat", line[2])
            if len(line) < 4:
                cf.set("Basis", "control", "")
            else:
                cf.set("Basis", "control", line[3])
                cf.set("Basis", "output", os.path.abspath(line[0]))
            cf.write(open(line[0] + ".conf", "w"))

            ## write in environment setting
            bash = """
#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -N python

#$ -pe mpi 1
#$ -j y
#$ -q ahost.q

source /software/intel/Compiler/11.1/038/bin/ifortvars.sh intel64
source /software/intel/Compiler/11.1/038/bin/iccvars.sh intel64
source /software/intel/mkl/10.2.0.013/tools/environment/mklvarsem64t.sh

source /usr/mpi/intel/mvapich-1.2.0/bin/mpivars.sh
#source /usr/mpi/gcc/mvapich-1.2.0/bin/mpivars.sh
#source /software/intel/impi/3.2.1.009/bin64/mpivars.sh

alias qstat='/home/script/zhang.sh'
#export PYTHONPATH=/software/python.2.7.3/lib/python2.7/site-packages:${PYTHONPATH}
export PATH=/software/python.2.7.3/bin:$PATH
export BLAS=~/src/BLAS/libfblas.a
export LAPACK=/homea3/linxq/src/lapack-3.4.1/libflapack.a
export ATLAS=/homea3/linxq/src/ATLAS/ATLAS/lib/libatlas.a
#export PATH=/install/rhels5.5/x86_64/Server/python-devel-2.4.3-27.el5.x86_64.rpm
#export LAPACK=/homea3/linxq/src/lapack-3.1.1/lapack_LINUX.a
#export PATH=/software/ge-6.1/bin/lx24-amd64:$PATH:/software/bin
source /software/ge-6.1/default/common/settings.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/ATLAS/ATLAS_x86/lib

###### my bashrc #######
export PATH=$PATH:/software/samtools-0.1.18
export PATH=/homea6/qinqian/opt/bin:$PATH
export PATH=/homea6/qinqian/opt:$PATH
export PATH=/homea6/qinqian/bin:$PATH
export PATH=/software/git-1.7.6/:$PATH
export PYTHONPATH=/homea6/qinqian/opt/lib/python2.7/site-packages/:${PYTHONPATH}
export PYTHONPATH=/homea6/qinqian/opt/lib/python3.3/site-packages/:${PYTHONPATH}
export PYTHONPATH=/homea6/qinqian/opt:${PYTHONPATH}
export PATH=/software/python.3.3.0/bin/:${PATH}

export PATH=$PATH:/software/python.2.7.3/bin:/homea3/linxq/soft:/homea2/jxfeng/soft/BEDTools-Version-2.16.2/bedtools/bin:/homea2/jxfeng/linxq/biosoft/bismark_v0.7.1

export PATH=$PATH:/software/R-2.14.1/bin/:/software/python.2.7.3/bin:/homea3/linxq/soft:/homea2/jxfeng/soft/BEDTools-Version-2.16.2/bedtools/bin:/homea2/jxfeng/linxq/biosoft/bismark_v0.7.1
export PATH=$PATH:/software/blat_34/bin/x86_64/:/software/bowtie_0.12.7:/software/cufflinks-1.2.0.Linux_x86_64:/software/bwa_0.6.1



#nproc=$NSLOTS
#cp $TMPDIR/machines nodefile

ChiLin2.py run -c %s.conf  --skip 7
            """ % line[0]
            open("chilin_"+line[0] + ".sh", "w").write(bash)

if __name__ == "__main__": main()
