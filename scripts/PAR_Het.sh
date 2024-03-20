bcftools view -g het -m 2 -M2  mm2_DV/Y_HAPLOID.beagle4.vcf.gz | bcftools +fill-tags -- -t all | bcftools query -f '%INFO/AC_Het\t%INFO/DR2\n' > het_cor.csv
