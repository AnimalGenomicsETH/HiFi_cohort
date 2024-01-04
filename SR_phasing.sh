
echo "sample chromosome N50"
for C in {1..29}
do 
  for S in $(cat /cluster/work/pausch/alex/Pop_HiFi/sample_IDs.txt)
  do
    tail -n +2 /cluster/work/pausch/naveen/phASER/phASER/CHR${C}/${S}.haplotypic_counts.txt | cut -f 1-3 > ${S}.${C}.SR.bed
    echo "${S} ${C} $(awk '{print NF"\t"$3-$2+1}' ${S}.${C}.SR.bed | calN50.js -L $(awk -v C=$C '$1==C {print $2}' /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai) - | awk '$1=="NL"&&$2==0 {print $3}')"
    #echo "${S} ${C} $(calN50.js -L $(awk -v C=$C '$1==C {print $2}' /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai) ${S}.${C}.SR.bed | awk '$1=="NN" {print $2}')"
  done
done
