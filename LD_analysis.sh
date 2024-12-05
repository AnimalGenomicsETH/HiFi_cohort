for i in {1..29}
do
  bcftools annotate --set-id "%VKX" -o ${i}.hex.vcf.gz -W=tbi filtered.${i}.vcf.gz
  zgrep "SVTYPE=" ${i}.hex.vcf.gz | cut -f 3 > ${i}.SV_list
  plink2 --vcf ${i}.hex.vcf.gz --out SVs_${i} --threads 2  --vcf-half-call missing --ld-snp-list ${i}.SV_list --r2-unphased zs --ld-window 999999999 --cow
done
