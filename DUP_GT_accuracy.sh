#!/bin/bash

#awk -v OFS='\t' 'NR<32 {print $1,1,$2}' /cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa.fai > regions.bed

#bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE' /cluster/work/pausch/HiFi_QTL/final_set/final_filtered.all.vcf.gz | grep "DUP" | cut -f -3 >> regions.bed

echo "sample read chromosome start end coverage count" | sed 's/ /\t/g' > duplications.csv

parallel --line-buffer -j $1 "samtools bedcov --reference /cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa -c regions.bed /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_HiFi/{}.mm2.cram | awk -v OFS='\t' -v S={} '{print S,\"HiFi\",\$0}'; samtools bedcov --reference /cluster/work/pausch/inputs/ref/BTA/UCD2.0/masked/GCA_002263795.4_ARS-UCD2.0_genomic.hard_masked_PAR_X.fa.gz -c regions.bed /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_SR/{}.bwa.cram | awk -v OFS='\t' -v S={} '{print S,\"SR\",\$0}'" ::: $(bcftools query -l /cluster/work/pausch/HiFi_QTL/final_set/final_filtered.all.vcf.gz) >> duplications.csv



#echo -e "chrom\tpos\t$(bcftools query -l final_filtered.all.vcf.gz | tr '\n' '\t')" | sed 's/\t$//' > duplication_GTs.csv
#bcftools query -i 'INFO/SVTYPE=="DUP"' -f '%CHROM\t%POS[\t%GT]' final_filtered.all.vcf.gz >> duplication_GTs.csv
