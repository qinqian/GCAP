#!/bin/bash

## extract 1bp base pairs from bam files
## bam to sam, sam awk
## sam is 1-based, should coordinate - 1, or simply use bedtools
## bed is 0-based

if [ $# -lt 3 ];then
    echo `basename $0` "can run hotspot pipeline v3 in batch"
    echo "Need 1 parameters now: <bam/bed_file> <tags_file> <sequence_type>"
    exit 1
fi

bam=$1
bed=$1
tags=$2
seqtype=$3

## get tags for __TAGS__ in hotspot
# -eq for integer, = for string
if [[ $seqtype = "se" ]]
then
    echo "$seqtype, single end"
    if [ ! -e $tags ]
    then
    samtools view $bam \
          | awk 'BEGIN{FS="\t"; OFS="\t"}{if ($2!=4) {len=length($10); s=$4-1; if(and(16,$2) == 16) s+=len-1; chr=$3; print chr,s,s+1}}' - \
          | sort-bed - \
          | starch - > $tags
    else
        echo "$tags Already exists!"
    fi
elif [[ $seqtype = "pe" ]]
then
    echo "$seqtype, pair end"
    if [ ! -e $tags ]
    then
       bedtools bamtobed -i  $bam | awk  'BEGIN{FS="\t"; OFS="\t"} {if ($6=="+") {print $1,int($2),int($2)+1 } else {print $1,int($3)-1,int($3)}}' - | sort-bed - | starch - > $tags
    else
       echo "$tags Already exists!"
    fi
elif [[ $seqtype = "bed" ]]
then
    echo "$seqtype, bed files"
    if [ ! -e $tags ]
    then
        ## calculate redundancy
        cat $bed | awk  'BEGIN{FS="\t"; OFS="\t"} {if ($6=="+") {print $1,int($2),int($2)+1 } else {print $1,int($3)-1,int($3)}}' - | sort-bed - | starch - > $tags
        unstarch $tags | tee ${tags}.tmp | wc -l > ${tags}.loc_count
        uniq ${tags}.tmp | wc -l > ${tags}.uniq_loc
        rm ${tags}.tmp
    else
        echo "$tags Already exists!"
    fi
fi