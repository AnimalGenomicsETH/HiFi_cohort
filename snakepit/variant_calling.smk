


rule deepvariant:
    input:
        expand(rules.samtools_merge.output,samples=samples)
    output:
        '{read_type}_DV/cohort.Unrevised.vcf.gz'
    shell:
        '''
        do
        '''

rule beagle4_impute:
    input:
        rule.deepvariant.output
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
        cram = phased + "{sample}/{sample}.mm2.haplotagged.cram",
        vcf_snv = phased + "{sample}/{sample}.mm2.DV.phased.vcf.gz",
        summary = multiext("{sample}/{sample}.mm2.phased",".summary.tsv","stats.tvt","blocks.tsv")
    threads: 16
    resources:
        mem_mb = 1500,
        walltime = '4h'
    shell:
        '''
        hiphase -t {threads} -r {input.reference} -b {input.cram} -p {output.cram} -c {input.vcf} -o {output.vcf} -s {wildcards.sample} --summary-file {output.summary[0]} --stats-file {output.summary[1]} --blocks-file {output.summary[2]}
        '''
