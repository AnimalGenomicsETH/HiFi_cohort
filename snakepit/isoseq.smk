

def get_SMRT_cells_per_sample(sample):
    return config['samples'][sample]

rule all:
    input:
        expand('isoseq/{sample}.gff',sample=config['samples'])

rule isoseq_refine:
    input: 
        bam = lambda wildcards: config['samples'][wildcards.sample][wildcards.SMRT]#'isoseq/{sample}/{SMRT}.bam',
    output:
        bam = 'isoseq/{sample}/{SMRT}.flnc.bam'
    params:
        primers = lambda wildcards: '\\n'.join(f'>{end}\\n{sequence}' for end,sequence in config['primers'][wildcards.sample].items())
    conda: 'pbccs'
    threads: 8
    resources:
        mem_mb = 1500
    shell:
        '''
        echo -e "{params.primers}" > $TMPDIR/primers.fa
        isoseq refine -j {threads} {input.bam} $TMPDIR/primers.fa {output.bam}
        '''

rule isoseq_cluster2:
     input:
         lambda wildcards: expand(rules.isoseq_refine.output['bam'],SMRT=get_SMRT_cells_per_sample(wildcards.sample),sample=wildcards.sample)
     output:
         bam = 'isoseq/{sample}/clustered.bam'
     conda: 'pbccs'
     threads: 8
     resources:
         mem_mb = 2500
     shell:
         '''
         isoseq cluster2 <(ls {input.bam}) {output.bam}
         '''

rule isoseq_pbmm2:
    input:
        bam = rules.isoseq_cluster2.output['bam'],
        reference = config['reference']
    output:
        bam = 'isoseq/{sample}/aligned.bam'
    conda: 'pbccs'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        '''
        pbmm2 align --preset ISOSEQ --sort {input.bam} {input.reference} {output.bam}
        '''

rule isoseq_collapse:
    input:
        mapped_bam = rules.isoseq_pbmm2.output['bam'],
        flnc_bam = lambda wildcards: expand(rules.isoseq_refine.output['bam'],SMRT=get_SMRT_cells_per_sample(wildcards.sample),sample=wildcards.sample)
    output:
        gff = 'isoseq/{sample}.gff'
    conda: 'pbccs'
    threads: 4
    resources:
        memb_mb = 5000
    shell:
        '''
        isoseq collapse --do-not-collapse-extra-5exons {input.mapped_bam} <(ls {input.flnc_bam}) {output.gff}
        '''
