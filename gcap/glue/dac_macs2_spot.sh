#!/bin/bash

#########################################################
#
# Compute SPOT (Signal Portion Of Tags) metric for MACS2
#
#########################################################

if [ $# -lt 2 ];then
    echo `basename $0` "calculate MACS2 SPOT score"
    echo "Need 1 parameters now: <tags bam/bed file> <peaks bed file> <out>"
    exit 1
fi

tags=$1
peaks=$2

thisscr="macs2_spot.sh"
echo
echo $thisscr

# Check tags for proper naming.
test=$(echo $tags | grep "\.bam$")
if [ ${#test} != 0 ]; then
    bam=T
else
    test=$(echo $tags | grep "\.bed$")
    if [ ${#test} != 0 ]; then
	    bam=F
    else
	    echo "$thisscr: $tags must end in .bam or .bed"
	    exit
    fi
fi

proj=`echo $tags | sed s/\.bam$// | sed s/\.bed$//`

echo ${proj}_tags.bed
if [ ! -e ${proj}_tags.bed ]
then
    if [ $bam == "T" ]; then
        bamToBed -i $tags \
            | awk 'BEGIN{OFS="\t"}{if($6 == "-") $2=$3-1; print $1, $2, $2+1}' \
            | sort-bed - > ${proj}_tags.bed
    else
        cat $tags \
            | awk 'BEGIN{OFS="\t"}{if($6 == "-") $2=$3-1; print $1, $2, $2+1}' \
            | sort-bed - > ${proj}_tags.bed
    fi
fi

ntag=`wc -l ${proj}_tags.bed | cut -d" " -f1`

out=$3

tih=$(cat ${proj}_tags.bed | bedops -e -1 - $peaks | wc -l)

spot=$(echo "scale=4; $tih/$ntag" | bc)

echo "tih = $tih"
echo "SPOT = $spot"

printf "%12s  %12s  %6s\n" "total tags" "MACS2 tags" "SPOT" > $out
printf "%12d  %12d  %.4f\n" $ntag $tih $spot >> $out
