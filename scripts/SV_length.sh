{ echo "sample length" ; for i in BSWCHEM1201*.sniffles.vcf.gz; do bcftools query -i 'ILEN>50&ILEN<50000' -f '%INFO/SVLEN\n' $i | awk -v S=${i%.sniffles.vcf.gz} '{print S,$1}'; done ; } > SV_lens.csv
