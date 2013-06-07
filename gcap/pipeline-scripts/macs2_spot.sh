#!/bin/bash

#########################################################
#
# Compute SPOT (Signal Portion Of Tags) metric for MACS2
# A little strange
#########################################################

if [ $# -lt 2 ];then
    echo `basename $0` "calculate MACS2 SPOT score"
    echo "Need 1 parameters now: <tags bed starch file> <peaks bed file>"
    exit 1
fi

tags=$1
peaks=$2

thisscr="macs2_spot.sh"
echo
echo $thisscr

ntag=`unstarch $tags | wc -l | cut -d" " -f2`

out=${peaks}.spot.out

tih=$(unstarch $tags | bedops -e -1 - $peaks | wc -l)

spot=$(echo "scale=4; $tih/$ntag" | bc)

echo "tih = $tih"
echo "SPOT = $spot"

printf "%12s  %12s  %6s\n" "total tags" "MACS2 tags" "SPOT" > $out
printf "%12d  %12d  %.4f\n" $ntag $tih $spot >> $out