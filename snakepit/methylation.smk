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
        cram = rules.hiphase.output['bam'],
        reference = config['reference']
    output:
        bed = "{sample}/{sample}.pbmm2.combined.bed",
        BigWig = "{sample}/{sample}.pbmm2.combined.bw"
    threads: 4
    resources:
        mem_mb = 1500,
        walltime = '4h'
    params:
        prefix = lambda wildcards, output: PurePath(output['bed']).with_suffix('').with_suffix(''),
        model = '/cluster/work/pausch/alex/software/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite'
    shell:
        '''
        aligned_bam_to_cpg_score --bam {input.cram} --ref {input.reference} --model {params.model} --output-prefix {params.prefix} --threads {threads}
        '''

rule methbat_profile:
    input:
        expand(rules.pb_CpG_tools.output['bed'],sample=samples)
    output:
        'methylation/merged.bed'
    shell:
        '''
        methbat profile --input-prefix {params.prefix} --output-region-profile {output}
        '''

rule methbat_build:
    input:
        collection = '',
        profiles = expand(rules.methbat_profile.output)
    output:
        ''
    shell:
        '''
        methbat build --input-collection {input.collection} --output-profile {output.profile}
        '''

rule methbat_compare:
    input:
        rules.methbat_build.output
    output:
        ''
    shell:
        '''
        methbat compare --input-profile {input} --output-comparison {output}
        '''
