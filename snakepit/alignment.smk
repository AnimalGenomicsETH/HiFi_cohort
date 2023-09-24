rule fibertools_predict_m6a:
    input:
        'alignments/uBAM/{sample}/{cell}.5mC.bam' #raw PacBio ubam with kinetics
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

rule minimap2_align:
    input:
        uBAM = lambda wildcards: rules.fibertools_predict_m6a.output[0] if wildcards.methylation == 'm6a' else 'alignments/uBAM/{sample}/{cell}.5mC.bam',
        reference = config['mm2_index']
    output:
        temp(multiext('alignments/mm2/{sample}/{cell}.{methylation}.bam','','.csi'))
    threads: lambda wildcards: 24 if wildcards.cell[:2] == 'm8' else 12
    resources:
        mem_mb = 3000,
        walltime = '24h',
        scratch = '50G'
    shell:
        '''
        samtools fastq -T rq,MM,ML --threads {threads} {input.uBAM} |\
        minimap2 -t {threads} -a map-hifi -R '@RG\\tPL:PacBio\\tID:{wildcards.cell}\\tSM:{wildcards.sample}' -y -Y {input.reference} - |\
        samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output[0]} --write-index
        '''

def gather_cells(sample):
    zipper = {'cell':[],'methylation':[]}
    for _,row in HiFi_dataframe[HiFi_dataframe['ID']==sample].iterrows():
        zipper['cell'].append(row["Metadata Context ID"])
        zipper['methylation'].append({"5mC" if row["Kinetic Data"} == "No" else "m6a"})
    return zipper

rule samtools_merge:
    input:
        bam = lambda wildcards: expand(rules.minimap2_align.output[0],zip,**gather_cells(wildcards.sample),allow_missing=True),
        csi = lambda wildcards: expand(rules.minimap2_align.output[1],zip,**gather_cells(wildcards.sample),allow_missing=True)
    output:
        multiext('alignments/{sample}.mm2.bam','','.csi')
    threads: 6
    resources:
        mem_mb = 5000
    shell:
        '''
        samtools merge -@ {threads} --write-index --reference {config[reference]} -c -o {output[0]} {input.bam}
        '''
