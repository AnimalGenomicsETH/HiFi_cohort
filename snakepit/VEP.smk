workflow._singularity_args += ' -B /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data -B /cluster/work/pausch/inputs '

rule VEP:
    input:
        vcf = 'final_set/final_filtered.all.vcf.gz',
        gff = '/cluster/work/pausch/HiFi_QTL/GCF_002263795.3_ARS-UCD2.0_genomic.gff',
        reference = '/cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa'
    output:
        vcf = multiext('/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/projects/eQTL_cohort/HiFi/refseq_vep.vcf.gz','','.tbi'),
        stats = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/projects/eQTL_cohort/HiFi/refseq.stats'
    container: '/cluster/work/pausch/alex/software/images/ensembl-vep_release_113.0.sif'
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '4h'
    shell:
        '''
        vep -i {input.vcf} -gff {input.gff} -fasta {input.reference} --hgvs --symbol -o {output.vcf[0]} --compress_output bgzip --stats_file {output.stats} --stats_text --vcf
        tabix -p vcf {output.vcf[0]}
        '''
