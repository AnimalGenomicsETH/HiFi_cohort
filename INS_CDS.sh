bcftools +split-vep -i 'INFO/SVTYPE="INS"' -f '%CHROM\t%POS\t%ID\t%IMPACT' VEP/refseq_vep.eGene.vcf.gz | awk -v OFS='\t' '{print $1,$2,$2+1,$3,$4}' > INS.bed

zcat VEP/ARS-UCD2.0_genomic.eGene.gtf.gz | awk -v OFS='\t' '$3=="CDS" {print $1,$4,$5,$10}' | sort -k1,1V -k2,2n | bedtools merge -d 0 -o distinct -c 4 > CDS.bed

bedtools intersect -wo -a CDS.bed -b INS.bed > SV_overlap.tsv
