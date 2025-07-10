# Snakemake to identify splice junctions (Regtools), map junctions to intron clusters (Leafcutter), and prepare for QTLtools 

# Get samples/chr
samples = [i.rstrip() for i in open('117_samples.txt')]

# Tools
regtools = '/path/to/regtools/build/regtools'
leafcutter_cluster = '/path/to/leafcutter/clustering/leafcutter_cluster_regtools.py'
clusters_to_genes = '/path/to/map_clusters_to_genes.R'
leafcutter_phenos = '/path/to/leafcutter/scripts/prepare_phenotype_table.py'
merge_sqtl = '/path/to/merge_sqtl.R'
combine_cov = '/path/to/combine_cov.R'

# Output
fold_out = '/cluster/work/pausch/HiFi_QTL/splicing/testis/removed/'

rule all:
	input:
		junc_files = expand(fold_out + 'junctions/{sample}.junc', sample=samples),
		all_jucs = expand(fold_out + 'all_juncs.txt'),
		numer_counts = expand(fold_out + 'clusters/hifi_testis_perind_numers.counts.gz'),
		gene_map = expand(fold_out + 'clusters/gene_map.tsv'),
		og_phenos = expand(fold_out + 'clusters/hifi_testis_perind.counts.gz.qqnorm_1'),
		all_qqnorm = expand(fold_out + 'clusters/hifi_testis_perind.counts.gz.qqnorm_all'),
		pheno_file = expand(fold_out + 'phenotypes/all.leafcutter.bed'),
		sqtl_cov = expand(fold_out + 'covariates/RNA_cov_testis_sqtl.txt')
			
rule regtools:
	input:
		bam = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/wasp_out/UCD2.0_masked/{sample}/{sample}.bam'
	output:
		juncs = fold_out + 'junctions/{sample}.junc'
	resources:
		mem_mb = 5000,
		walltime = '4h'
	shell:
		regtools + ' junctions extract -a 8 -m 50 -M 500000 -s RF {input.bam} -o {output.juncs}'

localrules: junc_ref
rule jun_ref:
	output:
		all_juncs = fold_out + 'all_juncs.txt'
	shell:
		'''
		for bam in `ls /path/to/junctions/*.junc`; do echo $bam >> {output.all_juncs}; done
		'''

rule leafcutter_cluster:
	input: 
		all_junc = fold_out + 'all_juncs.txt'
	output:
		fold_out + 'clusters/hifi_testis_perind_numers.counts.gz',
		counts = fold_out + 'clusters/hifi_testis_perind.counts.gz'
	params:
		prefix = 'hifi_testis'
	resources:
		mem_mb = 10000,
		time = '4h'
	envmodules:
		'gcc/8.2.0',
		'python/3.8.5',
		'htslib' 
	shell:
		'''
		cd /cluster/work/pausch/HiFi_QTL/splicing/testis/removed/clusters/
		python {leafcutter_cluster} -j {input.all_junc} -m 50 -o {params.prefix} -l 500000
		'''

rule clusters_to_genes:
	input:
		counts = fold_out + 'clusters/hifi_testis_perind_numers.counts.gz',
		exons = '/cluster/work/pausch/xena/transcriptome_inputs/ARS_UCD2.0_Y/leafcutter/ARS-UCD2.0.genes.exons.txt'
	output:
		fold_out + 'clusters/gene_map.tsv'
	params:
		prefix = 'gene_map.tsv',
		dir = fold_out + 'clusters/'
	resources:
		mem_mb = 3000,
		time = '4h'
	envmodules:
		'gcc/8.2.0',
		'gsl',
		'r/4.0.2'
	shell:
		'''
		Rscript {clusters_to_genes} {input.counts} {input.exons} {params.prefix} --output_dir {params.dir}
		'''

rule make_og_phenos:
	input:
		counts = fold_out + 'clusters/hifi_testis_perind.counts.gz'
	output:
		fold_out + 'clusters/hifi_testis_perind.counts.gz.qqnorm_1'
	resources:
		mem_mb = 10000,
		time = '4h'
	envmodules:
		'gcc/8.2.0',
		'python/3.8.5',
		'htslib' 
	shell:
		'''
		python {leafcutter_phenos} {input.counts} -p 10
		'''
		
localrules: combine_qqnorm
rule combine_qqnorm:
	input: 
		fold_out + 'clusters/hifi_testis_perind.counts.gz.qqnorm_1'
	output:
		all_qqnorm = fold_out + 'clusters/hifi_testis_perind.counts.gz.qqnorm_all'
	params:
		wkdir = fold_out + 'clusters/'
	shell:
		'''
		cd {params.wkdir}
		head -n 1 hifi_testis_perind.counts.gz.qqnorm_1 > header_all
		sed -i '1d' hifi_testis_perind.counts.gz.qqnorm_*
		cat header_all hifi_testis_perind.counts.gz.qqnorm_* > {output.all_qqnorm}
		'''
		
rule merge:
	input: 
		all_qqnorm = fold_out + 'clusters/hifi_testis_perind.counts.gz.qqnorm_all',
		gene_map = fold_out + 'clusters/gene_map.tsv'
	output:
		pheno_file = fold_out + 'phenotypes/all.leafcutter.bed',
		PCA_all = fold_out + 'covariates/splicing_pca_all.txt',
		PCA_cov = fold_out + 'covariates/pc_covariates.txt',
		var_explained = fold_out + 'covariates/var_explained.txt',
	resources:
		mem_mb = 20000,
		time = '4h'
	envmodules:
		'gcc/8.2.0',
		'r/4.2.2'
	shell:
		'''
		Rscript {merge_sqtl} \
		{input.all_qqnorm} \
		{input.gene_map} \
		{output.pheno_file} \
		{output.PCA_all} \
		{output.PCA_cov} \
		{output.var_explained}
		'''

rule covariates:
	input:
		age_rin = '/path/to/age_rin.txt',
		rna_pca = fold_out + 'covariates/pc_covariates.txt',
		geno_pcs = '/path/to/pruned_hifi.eigenvec'
	output:
		sqtl_cov = fold_out + 'covariates/RNA_cov_testis_sqtl.txt'
	resources:
		mem_mb = 10000,
		time = '4h'
	envmodules:
		'gcc/8.2.0',
		'gsl',
		'r/4.2.2'
	shell:
		'''
		Rscript {combine_cov} {input.age_rin} {input.rna_pca} {input.geno_pcs} {output.sqtl_cov}
		'''
