#!/bin/bash

### GET HETEROZYGOUS SNPS ###

for type in imputed
do
for chr in {1..29} X Y_HAPLOID Y_PAR
do
sbatch --time=04:00:00 --mem-per-cpu=5000 --wrap="bcftools view -m2 -M2 -v snps /cluster/work/pausch/HiFi_QTL/variants/Illumina/DeepVariant/$type/$chr.vcf.gz > /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/$chr.snps.vcf; bgzip /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/$chr.snps.vcf"
done
done

### CONCATENATE ### 

for type in filtered imputed
do
sbatch --time=04:00:00 --mem-per-cpu=5000 --wrap="bcftools concat /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/{1..29}.snps.vcf.gz /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/X.snps.vcf.gz /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/Y_PAR.snps.vcf.gz /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/Y_HAPLOID.snps.vcf.gz -O z -o /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/all.snps.vcf.gz"
done

### SPLIT VCF BY SAMPLES ###

data_folder=/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/wasp_out/UCD2.0_Y
for type in imputed
do
for sample_id in `ls ${data_folder}`
do
sbatch --time=04:00:00 --mem-per-cpu=5000 --wrap="vcftools --gzvcf /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/all.snps.vcf.gz --indv ${sample_id} --recode --stdout | bgzip -c > /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/heterozygous/${sample_id}.snps.vcf.gz"
done
done

### GET HETS ###

# ONLY FILTERED
for type in imputed
do
for sample_id in `ls ${data_folder}`
do
bcftools view -i 'GT="het"' /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/heterozygous/${sample_id}.snps.vcf.gz > /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/heterozygous/${sample_id}.snps.het.vcf
rm /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/heterozygous/${sample_id}.snps.vcf.gz
bgzip /cluster/work/pausch/xena/hifi_cohort/variants/snps/illumina/$type/heterozygous/${sample_id}.snps.het.vcf
done
done

