{if (/^[^@]/) {if ($2!="4"){ul[$2"\t"$3"\t"$4] += 1}}} END{ for (ull in ul) printf ull"\t"ul[ull]"\n"}
