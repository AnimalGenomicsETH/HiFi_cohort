
rule all:
    input:
        'coverage/regions.csv.gz'


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
        bam = lambda wildcards: 'alignments/{sample}.mm2.bam' if wildcards.mapper=='mm2' else 'alignments_old/{sample}.bwa.cram',
        windows = expand(rules.bedtool_makewindows.output,window=100000),
        reference = config['reference']
    output:
        'coverage/{sample}.{mapper}.{filtering,default|secondary|quality}.csv'
    params:
        flags = lambda wildcards: {'default':'0','secondary':'256','quality':'0 -Q 10'}[wildcards.filtering]
    threads: 1
    resources:
        mem_mb = 15000
    shell:
        '''
        samtools bedcov -c -g {params.flags} --reference {input.reference} {input.windows} {input.bam} |\
        awk '{{print "{wildcards.sample}","{wildcards.mapper}","{wildcards.filtering}",$1,$2,$4/($3-$2),$5}}' > {output}
        '''

rule bedtools_coverage:
    input:
        windows = expand(rules.bedtool_makewindows.output,window=100000),
        vcfs = '{mapper}_DV/all.beagle4.vcf.gz'#,region=(list(map(str,range(1,30))) + ['X','Y','MT']),allow_missing=True)
    output:
        'coverage/{sample}.{mapper}.variants.csv'
    resources:
        mem_mb = 15000
    shell:
        '''
        bedtools coverage -a {input.windows} -b <(bcftools query -f '%CHROM %POS\\n' {input.vcfs} | awk -v OFS='\\t' '{{print $1,$2,$2+1}}') |\
        awk '{{print "{wildcards.sample}","{wildcards.mapper}","variants",$1,$2,$4/($3-$2),$7}}' > {output}
        '''

rule gather_samples:
    input:
        expand(rules.samtools_bedcov.output,sample=config['samples'],mapper=('mm2','bwa'),filtering=('default','secondary','quality')),
         expand(rules.bedtools_coverage.output,sample=config['samples'],mapper=('mm2','bwa'))
    output:
        'coverage/regions.csv.gz'
    localrule: True
    shell:
        '''
        {{ echo "sample mapper filtering chromosome position coverage reads" ;  cat {input} ; }} |
        pigz -p 2 -c > {output}
        '''
