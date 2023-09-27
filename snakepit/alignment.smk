
rule minimap2_align:
    input:
        uBAM = lambda wildcards: rules.fibertools_predict_m6a.output[0] if wildcards.methylation == 'm6a' else 'alignments/uBAM/{sample}/{cell}.5mC.bam',
        reference = Path(config['reference']).with_suffix('.map_hifi.mmi')
    output:
        temp(multiext('alignments/mm2/{sample}/{cell}.{methylation}.bam','','.csi'))
    threads: lambda wildcards: 24 if wildcards.cell[:2] == 'm8' else 16
    resources:
        mem_mb = 5000,
        walltime = '24h',
        scratch = '50G'
    shell:
        '''
        samtools fastq -T rq,MM,ML --threads {threads} {input.uBAM} |\
        minimap2 -t {threads} -ax map-hifi -R '@RG\\tPL:PacBio\\tID:{wildcards.cell}\\tSM:{wildcards.sample}' -y -Y {input.reference} - |\
        samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output[0]} --write-index
        '''

rule pbmm2_align:
    input:
        uBAM = lambda wildcards: rules.fibertools_predict_m6a.output[0] if wildcards.methylation == 'm6a' else 'alignments/uBAM/{sample}/{cell}.5mC.bam',
        reference = Path(config['reference']).with_suffix('.CCS.mmi')
    output:
        temp(multiext('alignments/pbmm2/{sample}/{cell}.{methylation}.bam','','.csi'))
    threads: lambda wildcards, input: 24 if input.size_mb > 20e3 else 12
    resources:
        mem_mb = 3000,
        walltime = '4h',
        scratch = '50G'
    conda:
        'pbccs'
    shell:
        '''
        pbmm2 align {input.reference} {input.uBAM} {output[0]} --preset CCS --sort -j {threads} --sample {wildcards.sample} --sort-memory 3000M --bam-index BAI
        '''

def gather_cells(sample):
    zipper = {'cell':[],'methylation':[]}
    for _,row in HiFi_dataframe[HiFi_dataframe['Sample']==sample].iterrows():
        zipper['cell'].append(row["Cell"])
        zipper['methylation'].append("5mC" if row["Kinetics"] == "No" else "m6a")
    return zipper

rule samtools_merge:
    input:
        bam = lambda wildcards: expand(rules.minimap2_align.output[0],zip,**gather_cells(wildcards.sample),allow_missing=True) if wildcards.mapper == 'mm2' else expand(rules.pbmm2_align.output[0],zip,**gather_cells(wildcards.sample),allow_missing=True),
        csi = lambda wildcards: expand(rules.minimap2_align.output[1],zip,**gather_cells(wildcards.sample),allow_missing=True) if wildcards.mapper == 'mm2' else expand(rules.pbmm2_align.output[1],zip,**gather_cells(wildcards.sample),allow_missing=True)
    output:
        multiext('alignments/{sample}.{mapper}.bam','','.csi')
    threads: 6
    resources:
        mem_mb = 5000
    shell:
        '''
        samtools merge -@ {threads} --write-index --reference {config[reference]} -c -o {output[0]} {input.bam}
        '''
