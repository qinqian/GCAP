#!/bin/bash

## sampling by percentage
## percentage * total_reads = 5M

sam=$1
prog=$2
output=$3


total=$(grep -v @ $sam | wc -l | cut -f 2 -d" ")
percent=$(echo "scale=10; 5000000/$total" | bc)

echo $sam
echo $total
echo $percent
echo $prog
echo $output

java -Xmx2g -XX:ParallelGCThreads=4 -jar $prog I=$sam O=$output P=$percent VALIDATION_STRINGENCY=SILENT