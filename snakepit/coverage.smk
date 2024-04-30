rule bedtool_makewindows:
    input:
        fai = config['reference']+".fai"
    output:
        'coverage/windows.{window}.bed'
    localrule: True
    shell:
        '''
        bedtools makewindows -w {wildcards.window} -g <(cut -f -2 {input.fai}) > {output}
        '''

rule samtools_bedcov:
    input:
        bam = rules.samtools_merge.output[0], #lambda wildcards: f'/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_{"HiFi" if wildcards.mapper=="mm2" else "SR"}/{wildcards.sample}.{wildcards.mapper}.cram',
        windows = expand(rules.bedtool_makewindows.output,window=100000),
        reference = config['reference']
    output:
        'coverage/{sample}.{mapper}.{filtering,default|secondary|quality}.csv'
    params:
        flags = lambda wildcards: {'default':'','secondary':'-g 256','quality':'-d 2 -Q 5','quality_sup':'-d 2 -Q 5 -G 2048'}[wildcards.filtering]
    threads: 1
    resources:
        mem_mb = 15000
    shell:
        '''
        samtools bedcov -c {params.flags} --reference {input.reference} {input.windows} {input.bam} |\
        awk '{{print "{wildcards.sample}","{wildcards.mapper}","{wildcards.filtering}",$1,$2,$4/($3-$2),$5}}' > {output}
        '''

rule samtools_coverage:
    input:
        cram = lambda wildcards: f'/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_{"HiFi" if wildcards.mapper=="mm2" else "SR"}/{wildcards.sample}.{wildcards.mapper}.cram',
        reference = config['reference']
    output:
        'coverage/{sample}.{mapper}.coverage.csv'
    shell:
        '''
        samtools coverage --reference {input.reference} {input.cram}
        '''

rule bedtools_coverage:
    input:
        windows = expand(rules.bedtool_makewindows.output,window=100000),
        vcfs = expand('{mapper}_DV/{region}.{filtering}.vcf.gz',region=main_regions,allow_missing=True)
    output:
        'coverage/{sample}.{mapper}.{filtering,filtered|beagle4-filtered}.csv'
    resources:
        mem_mb = 15000
    shell:
        '''
        bcftools concat --naive-force {input.vcfs} |\
        bcftools query -f '%CHROM %POS\\n' |\
        awk -v OFS='\\t' '{{print $1,$2,$2+1}}' |\
        bedtools coverage -a {input.windows} -b /dev/stdin |\
        awk '{{print "{wildcards.sample}","{wildcards.mapper}","{wildcards.filtering}",$1,$2,$4/($3-$2),$7}}' > {output}
        '''

rule gather_samples:
    input:
        expand(rules.samtools_bedcov.output,sample=config['samples'],mapper=('mm2','bwa'),filtering=('default','secondary','quality')),
         expand(rules.bedtools_coverage.output,sample=config['samples'],mapper=('mm2','bwa'),filtering=('filtered','beagle4_filtered'))
    output:
        'coverage/regions.csv.gz'
    localrule: True
    shell:
        '''
        {{ echo "sample mapper filtering chromosome position coverage reads" ;  cat {input} ; }} |
        pigz -p 2 -c > {output}
        '''
