

rule samtools_fastq_HiFi:
    input:
        config['HiFi_reads'] + '/{sample}/{cell}.bam'
    output:
        temp('{sample}/FASTQ/{cell}.F1.fq.gz')
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        EXT=$(echo {input} | awk -F'.' '{{print $NF}}')
        
        if [[ "$EXT" != "gz" ]]
        then
          samtools fastq -0 {output} -T rq,MM,ML --threads {threads} {input}
        else
         ln -s {input} {output}
        fi
        '''

rule minimap2_align:
    input:
        reads = '{sample}/FASTQ/{cell}.{haplotype}.fq.gz'
    output:
        temp(multiext('alignments/mm2/{sample}/{cell}.{haplotype}.bam','','.csi'))
    threads: lambda wildcards: 24 if wildcards.haplotype != 'unknown' else 8
    resources:
        mem_mb = 3000,
        walltime = lambda wildcards: '24h' if wildcards.cell[:2] == 'm8' else '4h',
        scratch = '50G'
    params:
        mapping = ' '.join(config.get('mm2_params',[])),
        RG = "-R '@RG\\tPL:PacBio\\tID:{cell}\\tSM:{sample}'"
    shell:
        '''
        minimap2 -t {threads} -a {params.mapping} {params.RG} -y -Y {config[reference]} {input.reads} |\
        samtools sort - -m 3000M -@ 4 -T $TMPDIR -o {output[0]} --write-index
        '''

rule samtools_merge:
    input:
        bam = lambda wildcards: expand(rules.minimap2_align.output[0],cell=HiFi_cells[wildcards.sample],haplotype=get_sample_haplotypes(wildcards.sample),allow_missing=True),
        csi = lambda wildcards: expand(rules.minimap2_align.output[1],cell=HiFi_cells[wildcards.sample],haplotype=get_sample_haplotypes(wildcards.sample),allow_missing=True)
    output:
        multiext('alignments/{sample}.{mapper}.cram','','.crai')
    threads: 12
    resources:
        mem_mb = 5000
    params:
        fmt = ' '.join(f'--output-fmt-option {OPTION}' for OPTION in config.get('cram_options',[]))
    shell:
        '''
        samtools merge -@ {threads} --write-index {params.fmt} --reference {config[reference]} -c -o {output[0]} {input.bam}
        '''

rule cramino_stats:
    input:
        rules.samtools_merge.output[0]
    output:
        'alignments/{sample}.{mapper}.stats'
    threads: 4
    resources:
        mem_mb = 5000
    shell:
        '''
        cramino --reference {config[reference]} -t {threads} -s {wildcards.sample} {input} > {output}
        '''
