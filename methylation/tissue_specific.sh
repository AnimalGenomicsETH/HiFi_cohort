if [ $# -eq 0 ]
  then
    THRESHOLD=0.01
else
    THRESHOLD=$1
fi

echo "TPM threshold set to ${THRESHOLD}"

for T in epididymis testis vas
do
  awk -v OFS='\t' -v THRESHOLD=${THRESHOLD} 'NR>1 {S=0; for (i=6;i<=NF;i++) {S+=$i;N=(NF-6)}; if(S/N>THRESHOLD) {print $1,$2,$3,$4,S/N }}' ${T}.tsv | sort -k1,1V -k2,2n > ${T}.expressed.list
done

bedtools multiinter -i *.expressed.list | awk 'FNR==NR&&NR>1 {C[$1""$2""$3]=$4; next} { if (FNR>1) {print $1,$2,$3,C[$1""$2""$3],$5}}' testis.tsv - > tissue_specific_genes.tsv
