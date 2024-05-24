samtools faidx /cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa 18:26066814-26106814 > seq.fa
TRF seq.fa 2 7 7 80 10 50 2000 -h -ngs


for i in $(cat ../samples.txt); do trgt genotype --genome /cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa --repeats repeats.bed --reads /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_HiFi/$i.mm2.cram --output-prefix $i; done

for i in *gz; do tabix -p vcf $i; done

bcftools merge -m all -o merge.vcf.gz *.vcf.gz
