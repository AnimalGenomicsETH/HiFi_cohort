
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

rule winnowmap_align:
    input:
        uBAM = 'alignments/uBAM/{sample}/{cell}.5mC.bam',
        reference = config['reference'],
        repetitive_kmers = 'repetitive_k15.txt'
    output:
        temp(multiext('alignments/wm2/{sample}/{cell}.5mC.bam','','.csi'))
    threads: lambda wildcards: 24 if wildcards.cell[:2] == 'm8' else 16
    resources:
        mem_mb = 5000,
        walltime = '24h',
        scratch = '50G'
    shell:
        '''
        samtools fastq -T rq,MM,ML --threads {threads} {input.uBAM} |\
        winnowmap -t {threads} -ax map-pb -R '@RG\\tPL:PacBio\\tID:{wildcards.cell}\\tSM:{wildcards.sample}' -y -Y {input.reference} - |\
        samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output[0]} --write-index
        '''

rule pbmm2_align:
    input:
        uBAM = lambda wildcards: rules.fibertools_predict_m6a.output[0] if wildcards.methylation == 'm6a' else 'alignments/uBAM/{sample}/{cell}.5mC.bam',
        reference = Path(config['reference']).with_suffix('.CCS.mmi')
    output:
        temp(multiext('alignments/pbmm2/{sample}/{cell}.{methylation}.bam','','.bai'))
    params:
        index = lambda wildcards, output: PurePath(output[1]).suffix.upper()[1:]
    threads: lambda wildcards, input: 24 if input.size_mb > 20e3 else 12
    resources:
        mem_mb = 4000,
        walltime = '24h',
        scratch = '50G'
    conda:
        'pbccs'
    shell:
        '''
        pbmm2 align {input.reference} {input.uBAM} {output[0]} --preset CCS --sort -j {threads} --sample {wildcards.sample} --sort-memory 3000M --bam-index {params.index}
        '''

def gather_cells(sample):
    zipper = {'cell':[],'methylation':[]}
    for _,row in HiFi_dataframe[HiFi_dataframe['Sample']==sample].iterrows():
        zipper['cell'].append(row["Cell"])
        zipper['methylation'].append("5mC" if row["Kinetics"] == "No" else "m6a")
    return zipper

rule samtools_merge:
    input:
        bam = lambda wildcards: expand('alignments/{mapper}/{sample}/{cell}.{methylation}.bam',zip,**gather_cells(wildcards.sample),allow_missing=True),
        csi = lambda wildcards: expand('alignments/{mapper}/{sample}/{cell}.{methylation}.bam.csi',zip,**gather_cells(wildcards.sample),allow_missing=True),
        #bam = lambda wildcards: expand(rules.minimap2_align.output[0],zip,**gather_cells(wildcards.sample),allow_missing=True) if wildcards.mapper == 'mm2' else expand(rules.pbmm2_align.output[0],zip,**gather_cells(wildcards.sample),allow_missing=True),
        #csi = lambda wildcards: expand(rules.minimap2_align.output[1],zip,**gather_cells(wildcards.sample),allow_missing=True) if wildcards.mapper == 'mm2' else expand(rules.pbmm2_align.output[1],zip,**gather_cells(wildcards.sample),allow_missing=True)
    output:
        multiext('alignments/{sample}.{mapper,mm2|pbmm2|wm2}.bam','','.csi')
    threads: 6
    resources:
        mem_mb = 5000
    shell:
        '''
        samtools merge -@ {threads} --write-index --reference {config[reference]} -c -o {output[0]} {input.bam}
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

rule gather_alignment_stats:
    input:
        expand(rules.cramino_stats.output,sample=samples,allow_missing=True)
    output:
        'alignments/{mapper}.stats.csv'
    localrule: True
    shell:
        '''
        {{ echo "sample,cell,gigabases,reads,read N50,read length,QV" ; cat {input} ; }} > {output}
        '''


rule fastp_filter:
    input:
        expand(config['short_reads']+'{{sample}}_R{N}.fastq.gz',N=(1,2))
    output:
        fastq = temp(expand('fastq/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True))
    params:
        min_quality = config.get('fastp',{}).get('min_quality',15),
        unqualified = config.get('fastp',{}).get('unqualified',40),
        min_length  = config.get('fastp',{}).get('min_length',15)
    threads: 4
    resources:
        mem_mb = 2500
    shell:
        '''
        fastp -q {params.min_quality} -u {params.unqualified} -g --length_required {params.min_length} --thread {threads} -i {input[0]} -o {output.fastq[0]} -I {input[1]} -O {output.fastq[1]} --json /dev/null --html /dev/null
        '''

def generate_SR_aligner_command(wildcards,input):
    match wildcards.mapper:
        case 'strobe':
            return f'strobealign {input.reference} {input.fastq} --rg-id {wildcards.sample}'
        case 'bwa':
            return f'bwa-mem2 mem -R "@RG\\tID:{wildcards.sample}\\tCN:UNK\\tLB:{wildcards.sample}\\tPL:illumina\\tSM:{wildcards.sample}" -Y {input.reference} {input.fastq}'
        case _:
            raise('Unknown aligner')

rule short_read_align:
    input:
        fastq = expand(rules.fastp_filter.output,allow_missing=True),
        reference = config['reference']
    output:
        bam = multiext('alignments/{sample}.{mapper,strobe|bwa}.cram','','.crai'),
        dedup_stats = 'alignments/{sample}.{mapper}.dedup.stats'
    params:
        aligner_command = lambda wildcards, input: generate_SR_aligner_command(wildcards,input)
    threads: 16 #lambda wildcards: 24 if wildcards.mapper == 'bwa' else 12
    resources:
        mem_mb = 3000,
        scratch = '50g',
        walltime = '24h'
    shell:
        '''
        {params.aligner_command} -t {threads} |\
        samtools collate -u -O -@ {threads} - |\
        samtools fixmate -m -u -@ {threads} - - |\
        samtools sort -T $TMPDIR -u -@ {threads} |\
        samtools markdup -T $TMPDIR -S -@ {threads} --write-index -f {output.dedup_stats} --reference {input.reference} - {output.bam[0]}
        '''
