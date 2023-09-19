rule fibertools_predict_m6a:
    input:
        'alignments/uBAM/{sample}/{cell}.5mC.bam'
        #raw PacBio ubam with kinetics
    output:
        'alignments/uBAM/{sample}/{cell}.m6a.bam' #new ubam with methylation but without kinetics
    conda:
        'fiber'
    threads: 16
    resources:
        mem_mb = 1500,
        walltime = '24h'
    shell:
        '''
        ft m6a -t {threads} {input} {output}
        '''

rule hiphase:
    input:
        cram = "{sample}.mm2.cram",
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

rule pb_CpG_tools:
    input:
        cram = rules.hiphase.output['cram'],
        reference = ref
    output:
        bed = phased + "{sample}/{sample}.pbmm2.combined.bed",
        BigWig = phased + "{sample}/{sample}.pbmm2.combined.bw"
    threads: 12
    resources:
        mem_mb = 1500,
        walltime = '4h'
    params:
        prefix = lambda wildcards, output: PurePath(output['bed']).with_suffix('').with_suffix(''),
        model = "/cluster/work/pausch/audald/software/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite"
    shell:
        '''
        aligned_bam_to_cpg_score --bam {input.cram} --ref {input.reference} --model {params.model} --output-prefix {params.prefix} --threads {threads}
        '''

rule batMeth:
    input:
        expand(rules.pb_CpG_tools.output['bed'],sample=samples)
    output:
        'methylation/merged.bed'
    shell:
        '''
        methBat
        '''
