

rule deepvariant:
    input:
        expand(rules.samtools_merge.output,sample=samples),
        config = 'config/DV.yaml'
    output:
        '{read_type}_DV/all.Unrevised.vcf.gz'
    params:
        name = lambda wildcards, output: PurePath(output[0]).parent,
        model = lambda wildcards: 'WGS' if wildcards.read_type == 'SR' else 'PACBIO'
    localrule: True
    shell:
        '''
        snakemake -s /cluster/work/pausch/alex/BSW_analysis/snakepit/deepvariant.smk --configfile {input.config} \
        --config Run_name="{params.name}" bam_path="alignments/" bam_index=".csi" bam_name="{{sample}}.bam" model="{params.model}" \
        --profile "slurm/fullNT" --nolock
        '''

rule beagle4_impute:
    input:
        rules.deepvariant.output
    output:
        '{read_type}_DV/cohort.beagle4.vcf.gz'
    shell:
        '''
        do
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
