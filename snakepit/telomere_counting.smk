rule find_telomeres:
    input:
        bam = 'alignments/{sample}.mm2.cram'
    output:
        telomeres = 'telomeres/{sample}.telo'
    threads: 2
    resources:
        mem_mb = 5000
    shell:
        '''
        samtools fasta --threads {threads} --reference {config[reference]} -0 /dev/stdout {input.bam} |\
        seqtk telo -d 1000 -s 100 > {output}
        '''
