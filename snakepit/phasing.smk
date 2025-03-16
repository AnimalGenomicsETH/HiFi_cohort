rule split_chromosomes:
    input:
        rules.samtools_merge.output
    output:
        expand('phasing/splits/{{sample}}.{{mapper,mm2|pbmm2}}.{region}.{ext}',ext=('bam','bam.csi'),region=regions[:-1])
    params:
        regions = regions[:-1]
    threads: 2
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '2h'
    shell:
        '''
        for i in {params.regions}
        do
          samtools view -@ {threads} --write-index -o phasing/splits/{wildcards.sample}.{wildcards.mapper}.${{i}}.bam {input} $i
        done
        '''

rule hiphase_chromosome:
    input:
        #bam = expand(rules.samtools_merge.output,mapper='mm2',allow_missing=True),
        bam = expand('phasing/splits/{sample}.mm2.{{region}}.{ext}',ext=('bam','bam.csi'),sample=samples),#expand(rules.split_chromosomes.output,mapper='mm2',sample=samples,allow_missing=True),
        vcf = expand(rules.merge_QTL_variants.output,mapper='mm2',allow_missing=True),
        reference = config['reference']
    output:
        bam = expand('phasing/{sample}.{{region}}.mm2.phased.{ext}',ext=('bam','bam.bai'),sample=samples),
        vcf = expand("phasing/cohort.{{region}}.mm2.phased.{ext}",ext=('vcf.gz','vcf.gz.tbi')),
        summary = multiext("phasing/{sample}.{region}.mm2.phased",".summary.tsv",".stats.tsv",".blocks.tsv")
    params:
        IO_list = lambda wildcards, input, output: ' '.join(f'--bam {input.bam[I]} --output-bam {output.bam[I]} --sample-name {S}' for I,S in enumerate(samples))
    threads: 12
    resources:
        mem_mb_per_cpu = 5000,
        runtime = '24h'
    shell:
        '''
hiphase --threads {threads} --io-threads 6 --reference {input.reference} \
{params.IO_list} \
--vcf {input.vcf[0]} --output-vcf {output.vcf[0]} \
--global-realignment-cputime 300 \
--min-vcf-qual 20 \
--phase-singletons \
--summary-file {output.summary[0]} --stats-file {output.summary[1]} --blocks-file {output.summary[2]}
        '''

rule summarise_phasing_chromosome:
    input:
        expand('phasing/summary.{region}.txt',region=regions[:-1])
    output:
        'phasing/summary.csv'
    localrule: True
    shell:
        '''
        echo "sample_name chromosome num_variants num_heterozygous num_phased num_blocks basepairs_per_block_max block_ng50" > {output}
        for R in {input}
        do
          filename=$(basename -- "$R")
          filename="${{filename%.*}}"
          filename="${{filename##*.}}"
          awk -v I=$filename '$2==I {{print $1,$2,$3,$4,$5,$9,$19,$NF}}' $R  >> {output}
        done
        '''

rule concat_phases:
    input:
        expand(rules.hiphase_chromosome.output.vcf[0],region=regions[:-1])
    output:
        'phasing/cohort.chromosomes.phased.vcf.gz'
    localrule: True
    shell:
        '''
        bcftools concat --naive-force -o {output} {input}
        '''

rule samtools_merge_phased:
    input:
        expand('phasing/{{sample}}.{region}.mm2.phased.bam',region=regions[:-1])
    output:
        multiext('phasing/{sample}.mm2.phased.bam','','.csi')
    threads: 4
    resources:
        mem_mb_per_cpu = 2500
    shell:
        '''
        samtools merge -@ {threads} --write-index -o {output[0]} {input}
        '''

rule hiphase:
    input:
        bam = expand(rules.samtools_merge.output,mapper='mm2',allow_missing=True),
        small = rules.bcftools_concat_QTL.output, 
    output:
        bam = multiext("alignments/{sample}.mm2.haplotagged.bam",'','.bai'),
        small = multiext("variants/{sample}.mm2.DV.phased.vcf.gz",'','.tbi'),
        summary = multiext("variants/{sample}.mm2.phased",".summary.tsv",".stats.tsv",".blocks.tsv")
    threads: 8
    resources:
        mem_mb_per_cpu = 4000,
        runtime = '8h'
    shell:
        '''
        hiphase --threads {threads} --io-threads 6 --reference {config[reference]} \
        --bam {input.bam[0]} --output-bam {output.bam[0]} \
        --vcf {input.small[0]} --output-vcf {output.small[0]} \
        --sample-name {wildcards.sample} \
        --global-realignment-cputime 300 \
        --min-vcf-qual 20 \
        --phase-singletons \
        --summary-file {output.summary[0]} --stats-file {output.summary[1]} --blocks-file {output.summary[2]}
        '''