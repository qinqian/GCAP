intersectBed -b test.promoter -a ../test/result3/testid_treatment-both-passes/testid_treatment.hotspot.twopass.fdr0.01.merge.pks.bed  -wa | sort-bed - | uniq | awk '{i+=int($3-$2)} END {print i}'
intersectBed -b test.promoter -a ../test/result3/testid_treatment-both-passes/testid_treatment.hotspot.twopass.fdr0.01.merge.pks.bed  -wa | sort-bed - | uniq | wc -l
intersectBed -b test.promoter -a test.tagbed  -wa | sort-bed - | uniq | awk '{i+=int($3-$2)} END {print i}'

sort-bed ../test/result3/testid_treatment-both-passes/testid_treatment.hotspot.twopass.fdr0.01.merge.pks.bed | uniq | awk '{i+=int($3-$2)} END {print i}' 
wc -l ../test/result3/testid_treatment-both-passes/testid_treatment.hotspot.twopass.fdr0.01.merge.pks.bed

wc -l test.promoter

