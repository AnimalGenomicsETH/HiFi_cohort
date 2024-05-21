{ echo "chromosome,process,type,N_variants" ; bcftools query -l mm2_DV/1.Unrevised.vcf.gz ; } | tr '\n' ','
echo ""
for C in {1..29} X Y Y_PAR Y_HAPLOID
do
  bcftools query -i 'TYPE="snp"' -f '[%GT,]\n' mm2_DV/${C}.Unrevised.vcf.gz | awk -v C=${C} -F ',' '{++n;for (i=1;i<=120;++i) {a[i]+= $i=="./."?1:0}} END {printf C",Unrevised,SNP,"n","; for (i=1;i<=120;++i) {printf a[i]","}; printf "\n" }'
  bcftools query -i 'TYPE!="snp"' -f '[%GT,]\n' mm2_DV/${C}.Unrevised.vcf.gz | awk -v C=${C} -F ',' '{++n;for (i=1;i<=120;++i) {a[i]+= $i=="./."?1:0}} END {printf C",Unrevised,INDEL,"n","; for (i=1;i<=120;++i) {printf a[i]","}; printf "\n" }'
  for V in DEL INS DUP INV
  do
    bcftools query -r ${C} -i INFO/SVTYPE="\""${V}"\"" -f '[%GT,]\n' mm2_SVs/cohort.sniffles.forced.vcf.gz | awk -v C=${C} -v V=${V} -F ',' '{++n;for (i=1;i<=120;++i) {a[i]+= $i=="./."?1:0}} END {printf C",forced,"V","n","; for (i=1;i<=120;++i) {printf a[i]","}; printf "\n" }'
    bcftools query -r ${C} -i INFO/SVTYPE="\""${V}"\"" -f '[%GT,]\n' mm2_SVs/cohort.sniffles.denovo.vcf.gz | awk -v C=${C} -v V=${V} -F ',' '{++n;for (i=1;i<=120;++i) {a[i]+= $i=="./."?1:0}} END {printf C",denovo,"V","n","; for (i=1;i<=120;++i) {printf a[i]","}; printf "\n" }'
  done
done
