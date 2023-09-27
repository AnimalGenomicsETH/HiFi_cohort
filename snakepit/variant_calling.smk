
regions = list(map(str,range(1,30))) + ['X','Y','MT','unplaced']
rule deepvariant:
    input:
        expand(rules.samtools_merge.output,sample=samples,mapper='mm2'),
        config = 'config/DV.yaml'
    output:
        expand('{read_type}_DV/{region}.Unrevised.vcf.gz',region=regions,allow_missing=True)
    params:
        name = lambda wildcards, output: PurePath(output[0]).parent,
        model = lambda wildcards: 'WGS' if wildcards.read_type == 'SR' else 'PACBIO'
    localrule: True
    shell:
        '''
        snakemake -s /cluster/work/pausch/alex/BSW_analysis/snakepit/deepvariant.smk --configfile {input.config} \
        --config Run_name="{params.name}" model="{params.model}" \
        --profile "slurm/fullNT" --nolock
        '''

rule beagle4_impute:
    input:
        rules.deepvariant.output
    output:
        multiext('{read_type}_DV/{region}.beagle4.vcf.gz','','.tbi')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        name = lambda wildcards, output: PurePath(output[0]).name
    threads: 10
    resources:
        mem_mb = 4000,
        walltime = '24h'
    shell:
        '''
        java -jar -Xss25m -Xmx40G /cluster/work/pausch/alex/software/beagle.27Jan18.7e1.jar gl={input} nthreads={threads} out={params.prefix}
        mv {output[0]} $TMPDIR/{params.name}
        bcftools reheader -f {config[reference]}.fai -o {output[0]} $TMPDIR/{params.name}
        tabix -fp vcf {output[0]}
        '''

rule bcftools_concat:
    input:
        expand(rules.beagle4_impute.output[0],allow_missing=True)
    output:
        multiext('{read_type}_DV/all.beagle4.vcf.gz','','.tbi')
    shell:
        '''
        bcftools concat {input} > {output}
        '''

rule pbsv_discover:
    input:
        expand(rules.samtools_merge.output,mapper='pbmm2',allow_missing=True)
    output:
        'SVs/{sample}.svsig.gz'
    conda:
        'pbccs'
    threads: 1
    resources:
        mem_mb = 2500
    shell:
        '''
        pbsv discover --ccs {input[0]} {output}
        '''

rule pbsv_call:
    input:
        signatures = expand(rules.pbsv_discover.output,sample=samples)
    output:
        'SVs/cohort.pbsv.vcf'
    conda:
        'pbccs'
    threads: 8
    resources:
        mem_mb = 4000
    shell:
        '''
        pbsv call --hifi -j {threads} --max-ins-length 20k {config[reference]} {input.signatures} {output}
        '''

rule sniffles_call:
    input:
        rules.samtools_merge.output
    output:
        vcf = 'SVs/{sample}.vcf.gz',
	snf = 'SVs/{sample}.snf'
    shell:
        '''
        '''

rule sniffles_merge:
    input:
        expand(rules.sniffles_call.output['snf'],sample=samples)
    output:
        'SVs/cohort.SVs.vcf.gz'
    shell:
        '''
        sniffles merge ...
        '''

rule hiphase:
    input:
        bam = rules.samtools_merge.output,
        vcf_snv = 'vcf.gz',
        reference = config['reference']
    output:
        cram = "{sample}/{sample}.mm2.haplotagged.cram",
        vcf_snv = "{sample}/{sample}.mm2.DV.phased.vcf.gz",
        summary = multiext("{sample}/{sample}.mm2.phased",".summary.tsv","stats.tvt","blocks.tsv")
    threads: 16
    resources:
        mem_mb = 1500,
        walltime = '4h'
    shell:
        '''
        hiphase -t {threads} -r {input.reference} -b {input.cram} -p {output.cram} -c {input.vcf} -o {output.vcf} -s {wildcards.sample} --summary-file {output.summary[0]} --stats-file {output.summary[1]} --blocks-file {output.summary[2]}
        '''
