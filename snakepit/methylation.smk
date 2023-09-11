rule fibertools_predict_m6a:
    input:
        '' #raw PacBio ubam with kinetics
    output:
        '' #new ubam with methylation but without kinetics 
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
