rule fibertools_predict_m6a:
    input:
        'alignments/uBAM/{sample}/{cell}.5mC.bam'
        #raw PacBio ubam with kinetics
    output:
        'alignments/uBAM/{sample}/{cell}.m6a.bam' #new ubam with methylation but without kinetics
    threads: 16
    resources:
        mem_mb = 1500,
        walltime = '24h'
    shell:
        '''
        ft m6a -t {threads} {input} {output}
        '''

rule pb_CpG_tools:
    input:
        cram = rules.hiphase.output['cram'],
        reference = config['reference']
    output:
        bed = "{sample}/{sample}.pbmm2.combined.bed",
        BigWig = "{sample}/{sample}.pbmm2.combined.bw"
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
