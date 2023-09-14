rule minimap2_align:
    input:
	uBAM = lambda wildcards: rules.fibertools_predict_m6a.output[0],
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

def merge_cells(HiFi_dataframe,sample):
    for _,row in HiFi_dataframe[HiFi_dataframe['ID']==sample].iterrows():
        f'{row["Metadata Context ID"]}.{"5mC" if row["Kinetic Data"} == "No" else "m6a"}.bam'
	#iterate over pandas to get cells with kinetics
	return None

rule samtools_merge:
    input:
        bam = lambda wildcards: expand(rules.minimap2_align.output[0],cell=HiFi_cells[wildcards.sample],haplotype=get_sample_haplotypes(wildcards.sample),allow_missing=True),
        csi = lambda wildcards: expand(rules.minimap2_align.output[1],cell=HiFi_cells[wildcards.sample],haplotype=get_sample_haplotypes(wildcards.sample),allow_missing=True)
    output:
        multiext('alignments/{sample}.mm2.bam','','.csi')
    threads: 6
    resources:
        mem_mb = 5000
    params:
        fmt = ' '.join(f'--output-fmt-option {OPTION}' for OPTION in config.get('cram_options',[]))
    shell:
        '''
        samtools merge -@ {threads} --write-index {params.fmt} --reference {config[reference]} -c -o {output[0]} {input.bam}
        '''
