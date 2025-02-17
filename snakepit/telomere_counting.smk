rule all:
    input:
        expand('telomeres/{sample}.telo',sample=config['samples'])

match config.get('cohort','Mature'):
    case 'Mature':
        path = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_HiFi/'
    case 'Calf':
        path = 'alignments/'

rule find_telomeres:
    input:
        bam = path + '{sample}.mm2.cram'
    output:
        telomeres = 'telomeres/{sample}.telo'
    threads: 2
    resources:
        mem_mb = 5000,
        walltime = '1h'
    shell:
        '''
        samtools fasta --threads {threads} --reference {config[reference]} -0 /dev/stdout {input.bam} |\
        seqtk telo -d 1000 -s 100 - > {output}
        '''
