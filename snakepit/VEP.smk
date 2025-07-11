rule prepare_GTF:
    output:
        multiext('ARS-UCD2.0_genomic.gtf.gz','','.tbi')
    localrule: True
    shell:
        '''
        wget -qO - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.gtf.gz |
        zgrep -ve "Curated" -ve "miRNA" -ve "#" |
        sort -k1,1V -k4,4n -k5,5n -t$'\\t' |
        bgzip -c > {output[0]}

        tabix -p gff {output[0]}
        '''

rule filter_GTF:
    input:
        gtf = rules.prepare_GTF.output,
        eGenes = 'gene_expression/testis/UCD2.0_masked/117_samples/TPM_gene_testis_UNFILTERED.tsv'
    output:
        multiext('VEP/ARS-UCD2.0_genomic.eGene.gtf.gz','','.tbi')
    localrule: True
    shell:
        '''
        zgrep -wf <(awk 'NR>1 {{print $4}}' {input.eGenes}) {input.gtf[0]} |\
        bgzip -c > {output[0]}

        tabix -p gff {output[0]}
        '''

#TODO: wildcard for taking in the {bwa,mm2}_DV variants
rule VEP:
    input:
        vcf = 'mm2_DV/1.Unrevised.vcf.gz',
        gtf = lambda wildcards: rules.prepare_GTF.output if wildcards.set == 'genomic' else rules.filter_GTF.output,
        reference = config['reference']
    output:
        vcf = multiext('VEP/refseq_vep.{set}.vcf.gz','','.tbi'),
        stats = 'VEP/refseq.stats.{set}.txt'
    container: '/cluster/work/pausch/alex/software/images/ensembl-vep_release_113.0.sif'
    threads: 8
    resources:
        mem_mb = 1000,
        walltime = '4h'
    shell:
        '''
        vep --fork {threads} -i {input.vcf} -gtf {input.gtf[0]} -fasta {input.reference} --hgvs --symbol -o {output.vcf[0]} --compress_output bgzip --stats_file {output.stats} --stats_text --vcf
        tabix -p vcf {output.vcf[0]}
        '''

rule overlap_annotations:
    input:
        vep = rules.VEP.output['vcf'],
        gtf = lambda wildcards: rules.prepare_GTF.output if wildcards.set == 'genomic' else rules.filter_GTF.output
    output:
        multiext('VEP/{set}','.INS.bed','.CDS.bed','.overlap.bed')
    localrule: True
    shell:
        '''
        bcftools +split-vep -i 'INFO/SVTYPE="INS"' -f '%CHROM\\t%POS\\t%ID\\t%IMPACT' {input.vep[0]} |\
        awk -v OFS='\\t' '{{print $1,$2,$2+1,$3,$4}}' > {output[0]}

        zcat {input.gtf[0]} |\
        awk -v OFS='\\t' '$3=="CDS" {{print $1,$4,$5,$10}}' |\
        sort -k1,1V -k2,2n |\
        bedtools merge -d 0 -o distinct -c 4 > {output[1]}

        bedtools intersect -wo -a {output[1]} -b {output[0]} > {output[2]}
        '''
