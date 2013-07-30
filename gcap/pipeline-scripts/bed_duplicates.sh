#!/bin/bash


#########################################################
#
# Calculate naive redundancy for 5M bed reads files
#
#########################################################

if [ $# -lt 3 ];then
    echo `basename $0` "can run hotspot pipeline v3 in batch"
    echo "Need 1 parameters now: <bed_file> <output_file>"
    exit 1
fi

bed=$1
outputstarch=$2
tool=$3
metrics=$4
total=$(wc -l $bed | cut -f1 -d" ")
echo $total
unique=$(cut -f 1,2,3 $bed | sort-bed - | uniq | wc -l | cut -f1 -d" ")
echo $total, $unique
redundancy=$(echo "scale=5; 1-$unique/$total" | bc)
echo $redundancy > $metrics

if [[ $tool = "hotspot" ]]
then
    sort-bed $bed | starch - > $outputstarch  ## starch to support bed.starch for hotspot v4 run_badspot
fi