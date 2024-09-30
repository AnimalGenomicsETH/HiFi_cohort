def get_SMRT_cells_per_sample(sample):
    return config['samples'][sample]

rule all:
    input:
        expand('isoseq/{sample}.gff',sample=config['samples'])

rule isoseq_refine:
    input: 
        bam = lambda wildcards: config['samples'][wildcards.sample][wildcards.SMRT]#'isoseq/{sample}/{SMRT}.bam',
    output:
        bam = multiext('isoseq/{sample}/{SMRT}.flnc.bam','','.pbi'),
        summary = multiext('isoseq/{sample}/{SMRT}.flnc','.consensusreadset.xml','.filter_summary.report.json','report.csv')
    params:
        primers = lambda wildcards: '\\n'.join(f'>{end}\\n{sequence}' for end,sequence in config['primers'][wildcards.sample].items())
    conda: 'pbccs'
    threads: 4
    resources:
        mem_mb = 1000,
        walltime = '30m'
    shell:
        '''
        echo -e "{params.primers}" > $TMPDIR/primers.fa
        isoseq refine -j {threads} {input.bam} $TMPDIR/primers.fa {output.bam[0]}
        '''

rule isoseq_cluster2:
     input:
         bams = lambda wildcards: expand(rules.isoseq_refine.output['bam'][0],SMRT=get_SMRT_cells_per_sample(wildcards.sample),sample=wildcards.sample)
     output:
         fofn = temp('{sample}.fofn'),
         bam = 'isoseq/{sample}/clustered.bam'
     conda: 'pbccs'
     threads: 8
     resources:
         mem_mb = 1000,
         walltime = '1h'
     shell:
         '''
         ls {input.bams} > {output.fofn}
         isoseq cluster2 -j {threads} {output.fofn} {output.bam}
         '''

rule isoseq_pbmm2:
    input:
        bam = rules.isoseq_cluster2.output['bam'],
        reference = config['reference']
    output:
        bam = multiext('isoseq/{sample}/aligned.bam','','.bai')
    conda: 'pbccs'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        '''
        pbmm2 align -j {threads} --preset ISOSEQ --sort {input.bam} {input.reference} {output.bam[0]}
        '''

rule isoseq_collapse:
    input:
        mapped_bam = rules.isoseq_pbmm2.output['bam'][0],
        flnc_bams = lambda wildcards: expand(rules.isoseq_refine.output['bam'][0],SMRT=get_SMRT_cells_per_sample(wildcards.sample),sample=wildcards.sample)
    output:
        gff = 'isoseq/{sample}.gff'
    conda: 'pbccs'
    threads: 4
    resources:
        memb_mb = 5000
    shell:
        '''
        isoseq collapse -j {threads} --do-not-collapse-extra-5exons {input.mapped_bam} <(ls {input.flnc_bams}) {output.gff}
        '''
