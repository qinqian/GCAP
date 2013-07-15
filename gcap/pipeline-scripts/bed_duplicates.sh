#!/bin/bash


#########################################################
#
# Calculate naive redundancy for bed reads files
#
#########################################################

if [ $# -lt 2 ];then
    echo `basename $0` "can run hotspot pipeline v3 in batch"
    echo "Need 1 parameters now: <bed_file> <output_file>"
    exit 1
fi

bed=$1
output=$2
tool=$3
total=$(wc -l $bed)
unique=$(cut -f 1,2,3 $bed | sort-bed - | uniq | wc -l)
redundancy=$(echo "scale=5; 1-$uniq/$total" | bc)
echo $redundancy > ${bed}.dup_metrics

if [[ $tool = "hotspot" ]]
then
    starch $bed > $output
fi