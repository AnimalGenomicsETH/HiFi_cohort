for i in INS DEL DUP INV
do
  echo "$i $(zgrep "@" SV_sequences.trf.gz | grep -c $i)"
done
for C in ERV RTE
do
  for i in INS DEL DUP INV
    do 
      echo "$C $i $(zgrep $C SV_sequences.out.gz | grep $i | awk '($7-$6)>250 {print $5}' | sort -u | wc -l)"
  done
done

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID' mm2_SVs/cohort.sniffles.forced.vcf.gz | awk -v OFS='\t' '{print $1,$2,$3+1,$4}' | grep -P "^(\d|X|Y|MT)" > cohort.sniffles.forced.bed
bedtools intersect -wo -a GCA_002263795.4_ARS-UCD2.0_genomic.trf.bed -F 0.5 -b cohort.sniffles.forced.bed | cut -f 7 | sort -u | grep -oP "(DEL|INS|INV|DUP)" | wc -l
