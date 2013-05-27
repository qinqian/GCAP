#!/bin/bash

#########################################################
#
# Compute SPOT (Signal Portion Of Tags) metric for MACS2
#
#########################################################

if [ $# -lt 3 ];then
    echo `basename $0` "can run hotspot pipeline v3 in batch"
    echo "Need 1 parameters now: <bam_file> <tags_file> <sequence_type>"
    exit 1
fi

## Not finish

tags=_TAGS_
outdir=_OUTDIR_

thisscr="macs2_spot.sh"
echo
echo $thisscr

proj=`basename $tags | sed s/\.bam$// | sed s/\.bed.starch$//`
tagb=$outdir/$proj.bed.starch
ntag=`cut -d" " -f2 $outdir/$proj-pass1/*.stdout`
hot=$outdir/$proj-both-passes/$proj.hotspot.twopass.zscore.wig
out=$outdir/$proj.spot.out
tih=$(unstarch $tagb | bedops --header -e -1 - $hot | wc -l)
spot=$(echo "scale=4; $tih/$ntag" | bc)
echo "tih = $tih"
echo "SPOT = $spot"
printf "%12s  %12s  %6s\n" "total tags" "hotspot tags" "SPOT" > $out
printf "%12d  %12d  %.4f\n" $ntag $tih $spot >> $out