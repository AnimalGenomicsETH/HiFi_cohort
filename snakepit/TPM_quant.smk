import os

# Kallisto TPM quantification and processing

# Get samples 
sample_file = '117_samples.txt'
sample_file = open(sample_file, 'r')
samples = [value.strip().split(' ')[0] for value in sample_file.readlines()]

# Tools
kallisto = '/path/to/kallisto/build/src/kallisto'

# Input
transcript_input = '/cluster/work/pausch/xena/transcriptome_inputs/masked_X/'
fastq = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/fast_p/'
tpm = '/cluster/work/pausch/HiFi_QTL/gene_expression/testis/UCD2.0_masked/'
kallisto_dir = '/cluster/work/pausch/HiFi_QTL/gene_expression/testis/UCD2.0_masked/kallisto/'

rule all:
	input:
		kallisto_quant = expand(kallisto_dir + '{sample}/abundance.tsv', sample=samples),
		TPM_mat = tpm + 'TPM_transcript_testis_FILTERED_NORM.tsv',
		final_cov = tpm + 'covariates/RNA_cov_testis.txt'
		
rule kallisto_quant:
	input:
		index = transcript_input + 'kallisto/ARS-UCD2.0_Xmasked_.refseq.cdna.idx',
		R1 = fastq + '{sample}_R1.fastq.gz',
		R2 = fastq + '{sample}_R2.fastq.gz'
	output:
		kallisto_dir + '{sample}/abundance.tsv'
	params:
		outdir = kallisto_dir + '{sample}'
	threads: 4
	resources:
		mem_mb = 2000,
		walltime = '4h'
	shell:
		kallisto + ' quant -i {input.index} -o {params.outdir} {input.R1} {input.R2} --rf-stranded -t {threads}'

rule make_matrix:
	input:
		kallisto_table = transcript_input + 'kallisto/tx2gene_table_cdna',
		gene_locations = '/cluster/work/pausch/xena/transcriptome_inputs/masked_X/gene_locos.tsv'
	output:
		TPM_unfiltered_transcript = tpm + 'TPM_transcript_testis_UNFILTERED.tsv',
		counts_unfiltered_transcript = tpm + 'counts_transcript_testis_UNFILTERED.tsv',
		TPM_unfiltered_gene = tpm + 'TPM_gene_testis_UNFILTERED.tsv',
		counts_unfiltered_gene = tpm + 'counts_gene_testis_UNFILTERED.tsv',
		TPM_filtered_transcript = tpm + 'TPM_transcript_testis_FILTERED_UNNORM.tsv',
		counts_filtered_transcript = tpm + 'counts_transcript_testis_FILTERED_UNNORM.tsv',
		TPM_normalized_transcript = tpm + 'TPM_transcript_testis_FILTERED_NORM.tsv',
		TPM_filtered_gene = tpm + 'TPM_gene_testis_FILTERED_UNNORM.tsv',
		counts_filtered_gene = tpm + 'counts_gene_testis_FILTERED_UNNORM.tsv',
		TPM_normalized_gene = tpm + 'TPM_gene_testis_FILTERED_NORM.tsv',
		PCs_all_transcript = tpm + 'covariates/PCs_all_transcript.txt',
		PCs_cov_transcript = tpm + 'covariates/PCs_cov_transcript.txt',
		var_explained_transcript = tpm + 'covariates/var_explained_transcript.txt',
		PCs_all_gene = tpm + 'covariates/PCs_all_gene.txt',
		PCs_cov_gene = tpm + 'covariates/PCs_cov_gene.txt',
		var_explained_gene = tpm + 'covariates/var_explained_gene.txt'
	params:
	    kallisto_dir = directory(tpm + 'kallisto'),
	    TPM_filter = '0.1',
	    read_filter = '6',
	    samp_filter = '0.2'
	envmodules:
		'stack/2024-06',
		'r/4.3.2'
	resources:
		mem_mb = 10000,
		walltime = '4h'
	shell:
		'''
		Rscript make_TPM.R \
		{input.kallisto_table} \
		{input.gene_locations} \
		{params.kallisto_dir} \
		{output.TPM_unfiltered_transcript} \
		{output.counts_unfiltered_transcript} \
		{output.TPM_unfiltered_gene} \
		{output.counts_unfiltered_gene} \
		{output.TPM_filtered_transcript} \
		{output.counts_filtered_transcript} \
		{output.TPM_filtered_gene} \
		{output.counts_filtered_gene} \
		{output.TPM_normalized_transcript} \
		{output.TPM_normalized_gene} \
		{params.TPM_filter} \
		{params.read_filter} \
		{params.samp_filter} \
		{output.PCs_all_transcript} \
		{output.PCs_cov_transcript} \
		{output.var_explained_transcript} \
		{output.PCs_all_gene} \
		{output.PCs_cov_gene} \
		{output.var_explained_gene}
		'''
	
rule final_cov:
	input:
		age_rin = '/cluster/work/pausch/xena/hifi_cohort/RNA/TPM_quantification/testis/covariates/age_rin.txt',
		rna_pcs = '/cluster/work/pausch/HiFi_QTL/gene_expression/testis/UCD2.0_masked/covariates/PCs_cov_gene.txt',
		geno_pcs = '/cluster/work/pausch/xena/hifi_cohort/variants/merged/pruned_hifi.eigenvec'
	output:
		final_cov = tpm + 'covariates/RNA_cov_testis.txt'
	envmodules:
		'stack/2024-06',
		'r/4.3.2'
	resources:
		mem_mb = 1000,
		walltime = '4h'
	shell:
		'''
		Rscript combine_cov.R {input.age_rin} {input.rna_pcs} {input.geno_pcs} {output.final_cov}
		'''