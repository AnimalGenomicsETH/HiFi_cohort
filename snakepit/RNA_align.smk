import os

#tissues = ['testis']
ref = ['UCD2.0_Y']

sample_file = "117_samples.txt"
sample_file = open(sample_file, 'r')
samples = [value.strip().split(' ')[0] for value in sample_file.readlines()] #only the first column is taken which has bam id

### TOOLS
STAR = '/path/to/STAR/source/STAR'
FASTP = '/path/to/fastp'
multiqc = '/path/to/multiqc'
FASTQ_SPLITTER = '/path/to/gdc-fastq-splitter'
SAMTOOLS = '/path/to/samtools/samtools'
PICARD_TOOLS = 'java -jar /cluster/apps/gcc-4.8.5/picard-2.25.7-krxq7vc5hgwrxwunjcwuefbjuzuvwawq/bin/picard.jar'
SAMBAMBA = '/path/to/sambamba_v0.6.6'

### INPUT AND OUTPUT
raw_data = '/path/to/raw_reads/'
vcf_path = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/het_vcfs/'
multiqc = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/multiqc/'
fast_p_out = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/fast_p/'
split_fastq = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/split_fastq/'
alignment = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/alignment/'
sorted_alignment = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/sorted_alignment/'
dedup_alignment = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/dedup_alignment/'
wasp_out = '/cluster/work/pausch/HiFi_QTL/RNA_alignment/testis/wasp_out/'

rule all:
	input:
		multiqc = expand(multiqc + "multiqc_report.html"),
		index= expand(dedup_alignment + "{ref}/{sample}/{sample}.bai", ref=ref, sample=samples),
		wasp = expand(wasp_out + "{ref}/{sample}/{sample}.bam", ref=ref, sample=samples),
        stats_1=expand(wasp_out + "{ref}/{sample}/{sample}.stats", ref=ref, sample=samples)

rule fastp:
	input:
		R1 = raw_data + "{sample}_R1.fastq.gz",
		R2 = raw_data + "{sample}_R2.fastq.gz"
	output:
		R1_O = temp(fast_p_out + "{sample}/{sample}_R1.fastq.gz"),
		R2_O = temp(fast_p_out + "{sample}/{sample}_R2.fastq.gz"),
		R_HTML = fast_p_out + "{sample}/{sample}_fastp.html",
		R_JSON = fast_p_out + "{sample}/{sample}_fastp.json"
	resources:
		mem_mb = 8000,
		time = "01:00:00"
	threads: 
		1
	shell:
		FASTP + " -i {input.R1} -o {output.R1_O} -I {input.R2} -O {output.R2_O} -h {output.R_HTML} -j {output.R_JSON} --trim_poly_g --trim_poly_x >/dev/null"

rule multiqc:
	input:
		fastp_json = expand(fast_p_out + "{sample}/{sample}_fastp.json",ref=ref,sample=samples)
	output:
		multiqc + "multiqc_report.html"
	resources:
		mem_mb = 1000,
		time = "01:00:00"
	threads: 
		1
	shell:
		multiqc + " -k json" + fast_p_out + "* -o" + fast_p_out + "multiqc"

checkpoint split_fastq:
	input:
		R1 = fast_p_out + "{sample}/{sample}_R1.fastq.gz",
		R2 = fast_p_out + "{sample}/{sample}_R2.fastq.gz"
	output:
		split_fastq = temp(directory(split_fastq + "{sample}"))
	resources:
		mem_mb = 5000,
		time = "24:00:00"
	threads:
		6
	params:
		prefix = split_fastq + "{sample}/{sample}_",
		folder = split_fastq + "{sample}"
	envmodules:
		'gcc/6.3.0',
		'python/3.7.4'
	shell:
		"mkdir {params.folder} \n" +
		FASTQ_SPLITTER + " -o {params.prefix} {input.R1} {input.R2}"

rule star:
	input:
		R1 = split_fastq + "{sample}/{sample}_{flowcell}_{lane}_R1.fq.gz",
		R2 = split_fastq + "{sample}/{sample}_{flowcell}_{lane}_R2.fq.gz",
		gDir = "/cluster/work/pausch/xena/transcriptome_inputs/ARS_UCD2.0_Y/STAR",
		het_vcf = vcf_path +  "{sample}.snps.het.vcf.gz",
		sjdb = "/cluster/work/pausch/xena/transcriptome_inputs/ARS_UCD2.0_Y/ARS-UCD2.0_genomic.gtf"
	output:
		bam = alignment + "{ref}/{sample}/{sample}_{flowcell}_{lane}.bam", 
		bai = alignment + "{ref}/{sample}/{sample}_{flowcell}_{lane}.bam.bai",
	resources:
		mem_mb = 7000,
		time = "04:00:00"
	threads:
		15
	params:
		rg = "SO=coordinate RGID={flowcell}:{lane} RGLB={sample}.0 RGPL=ILLUMINA RGPU=hiseq RGSM={sample} RGCN=FGCZ ",
	shell:
		"""
		cd $TMPDIR
		module load jdk 
		       
		/path/to/STAR/source/STAR --runThreadN 15 --twopassMode Basic --genomeDir {input.gDir} --sjdbGTFfile {input.sjdb} --sjdbOverhang 100 --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat --outSAMmapqUnique 60 --waspOutputMode SAMtag --varVCFfile <(zcat {input.het_vcf} ) --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random
				
		##coordinate sorting and add read groups
		java -jar /cluster/apps/gcc-4.8.5/picard-2.25.7-krxq7vc5hgwrxwunjcwuefbjuzuvwawq/bin/picard.jar \
		AddOrReplaceReadGroups I=Aligned.sortedByCoord.out.bam  O=Aligned_new_sorted.bam {params.rg} 

		## index 
		/path/to/samtools/samtools index Aligned_new_sorted.bam

		## copy files
		cp $TMPDIR/Aligned_new_sorted.bam {output.bam}
		cp $TMPDIR/Aligned_new_sorted.bam.bai {output.bai}
		"""    

def list_files(wildcards):
	checkpoint_output = checkpoints.split_fastq.get(sample = wildcards.sample).output.split_fastq
	all_wildcards = glob_wildcards(os.path.join(checkpoint_output, "{sample}_{flowcell}_{lane}_R1.fq.gz"))
	all_files = []
	for sample, flowcell, lane in zip(all_wildcards.sample, all_wildcards.flowcell, all_wildcards.lane):
		all_files.append(f"{alignment}" + "{ref}" + f"/{sample}/{sample}_{flowcell}_{lane}.bam")
	return(all_files)

rule sambamba_merge:
	input: list_files
	output: temp(sorted_alignment + "{ref}/{sample}/{sample}.bam")
	resources:
		mem_mb = 5000,
		time = "04:00:00"
	threads:
		6
	run:
		if len(input) == 1:
			shell("mv {input} {output} \n" + "mv {input}.bai {output}.bai")
		else:
			shell(SAMBAMBA + " merge -t {threads} {output} {input}")

rule sambamba_flagstat_1:
	input:
		BAM = sorted_alignment + "{ref}/{sample}/{sample}.bam"
	output:
		stats = sorted_alignment + "{ref}/{sample}/{sample}.stats"
	resources:
		mem_mb = 2000,
		time = "04:00:00"
	threads:
		10
	shell:
		SAMBAMBA + " flagstat -t {threads} {input.BAM} > {output.stats}"

rule mark_duplicates:
	input:
		sorted_alignment + "{ref}/{sample}/{sample}.bam"
	output:
		bam = dedup_alignment + "{ref}/{sample}/{sample}.bam",
		metrics = dedup_alignment + "{ref}/{sample}/{sample}.metrics.txt"
	resources:
		mem_mb = 8000,
		time = "04:00:00"
	threads:
		4
	envmodules:
		'gcc/6.3.0',
		'picard/2.25.7',
		'jdk'
	shell:
		PICARD_TOOLS + " MarkDuplicates I={input} O={output.bam} M={output.metrics}"

rule build_index:
	input:
		dedup_alignment + "{ref}/{sample}/{sample}.bam"
	output:
		dedup_alignment + "{ref}/{sample}/{sample}.bai"
	resources:
		mem_mb = 4000,
		time = "04:00:00"
	threads:
		6
	envmodules:
		'gcc/6.3.0',
		'picard/2.25.7',
		'jdk'
	shell:
		PICARD_TOOLS + " BuildBamIndex I={input} O={output}"

rule wasp_unique:
	input:
		BAM = dedup_alignment + "{ref}/{sample}/{sample}.bam",
		BAI = dedup_alignment + "{ref}/{sample}/{sample}.bai"
	output:
		BAM = wasp_out + "{ref}/{sample}/{sample}.bam"
	resources:
		mem_mb = 5000,
		time = "04:00:00"
	threads:
		1
	shell:
		"""
		module load python
		
		python  /path/to/Filter_wasp_unique.py {input.BAM} {output.BAM}
		"""
        
rule build_index2:
	input:
		wasp_out + "{ref}/{sample}/{sample}.bam"
	output:
		wasp_out + "{ref}/{sample}/{sample}.bam.bai"
	resources:
		mem_mb = 5000,
		time = "04:00:00"
	threads:
		1
	envmodules:
		'gcc/6.3.0',
		'picard/2.25.7',
		'jdk'
	shell:
		PICARD_TOOLS + " BuildBamIndex I={input} O={output}"

rule sambamba_flagstat_2:
	input:
		BAM = wasp_out + "{ref}/{sample}/{sample}.bam"
	output:
		stats = wasp_out + "{ref}/{sample}/{sample}.stats"
	resources:
		mem_mb = 2000,
		time = "04:00:00"
	threads:
		10
	shell:
		SAMBAMBA + " flagstat -t {threads} {input.BAM} > {output.stats}"
