BEGIN {FS="\t"; OFS="\t"}; {if (/^[^@]/) {total+=1; if ($2!="4" && $3!="chrM" && $3!="*" && $3!="chrY" && $3!="chrX" && $5>=1 && /XT:A:U/){ul+=1}}} END{ print ul"\t"total }