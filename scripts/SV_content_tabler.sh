{ echo "ID bases repeat_element repeat_class" ; awk '$5~/Sniffles2/ {print $5,($7-$6),$10,$11}' {1..29}.fa.out ; } > SV_repeats.csv

{ echo "ID motif number similarity score" ; awk '{ if($1~/@/) {C=substr($1,2)} else { if($1~/[[:digit:]+]/&&$3>=6&&$4>5) {print C,$14,$4,$6,$8 } } }' {1..29}.fa.trf | sort -k1,1 -k3,3n | awk '$5>a[$1] {a[$1]=$5;b[$1]=$0} END {for (k in b) {print b[k]}}' ; } > SV_VNTRs.csv
