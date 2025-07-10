/cluster/work/pausch/xena/programs/kallisto/build/src/kallisto  index -i /cluster/work/pausch/xena/transcriptome_inputs/ARS_UCD2.0_Y/kallisto/Bos_taurus.ARS-UCD2.0_Y.refseq.cdna.idx /cluster/work/pausch/xena/transcriptome_inputs/ARS_UCD2.0_Y/cow.1.rna.gz

echo "transcript_id gene_id gene_name" > tx2gene_table_cdna

zcat /cluster/work/pausch/xena/transcriptome_inputs/ARS_UCD2.0_Y/cow.1.rna.gz | grep ">" | awk -F")," '{print $(NF-2), $(NF-1)}' | awk -F"[(, ]" '{print $1, $NF, $NF}' | sed 's/>//g' >> tx2gene_table_cdna

