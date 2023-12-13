from pathlib import PurePath

wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    tissue = r'Testis',
    chunk = r'\d+',
    chrom = r'\d+',
    MAF = r'\d+',
    vcf = r'(eQTL|gwas)/\S+'

rule normalise_vcf:
    input:
        lambda wildcards: config['variants'][wildcards.variants]
    output:
        'QTL/{variants}.vcf.gz'
    threads: 4
    resources:
        mem_mb = 2500,
        disk_scratch = 50
    shell:
        '''
        bcftools norm --threads {threads} -f {config[reference]} -m -any {input} -Ou | \
        bcftools sort -T $TMPDIR -Ou - | \
        bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' -o {output} -
        tabix -fp vcf {output}
        '''

rule exclude_MAF:
    input:
        lambda wildcards: config['variants'][wildcards.variants]
    output:
        'QTL/{variants}.exclude_sites.{MAF}.txt'
    shell:
        '''
        bcftools view --threads 2 -Q 0.{wildcards.MAF}:minor -Ou {input} |\
        bcftools query -f '%ID\n'  - > {output}
        '''

def get_pass(_pass,input):
    if _pass == 'permutations':
        return f'--permute {config["permutations"]}'
    elif _pass == 'conditionals':
        return f'--mapping {input.mapping}'
    elif _pass == 'nominals':
        return f'--nominal {config.get("nominal",0.05)}'

rule qtltools_parallel:
    input:
        vcf = lambda wildcards: config['variants'][wildcards.variants],
        exclude = rules.exclude_MAF.output,
        bed = lambda wildcards: config['mol_QTLs'][wildcards.QTL][wildcards.tissue],
        cov = lambda wildcards: config['covariates'][wildcards.QTL][wildcards.tissue],
        mapping = lambda wildcards: 'QTL/{QTL}/{tissue}_{variants}/permutations_all.{MAF}.thresholds.txt' if wildcards._pass == 'conditionals' else []
    output:
        merged = temp('QTL/{QTL}/{tissue}_{variants}/{_pass}.{chunk}.{MAF}.txt.gz')
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input),
        grp = lambda wildcards: '--grp-best' if wildcards.QTL == 'sQTL' else ''
    threads: 1
    resources:
        mem_mb = 2500,
    shell:
        '''
        QTLtools cis --vcf {input.vcf} --exclude-sites {input.exclude} --bed {input.bed} --cov {input.cov} --std-err {params._pass} {params.grp} --window {config[window]} --normal --chunk {wildcards.chunk} {config[chunks]} --silent --log /dev/stderr --out /dev/stdout | pigz -p 2 -c > {output}
        '''

rule qtltools_gather:
    input:
        expand(rules.qtltools_parallel.output,chunk=range(0,config['chunks']+1),allow_missing=True)
    output:
        'QTL/{QTL}/{tissue}_{variants}/{_pass}.{MAF}.txt.gz'
    params:
        sort_key = lambda wildcards: '-k9,9n -k10,10n' if wildcards.QTL == 'eQTL' else '-k11,11n -k12,12n'
    localrule: True
    shell:
        '''
        LC_ALL=C; pigz -p 2 -dc {input} | sort --parallel=2 {params.sort_key} | pigz -p 2 -c > {output}
        '''

rule qtltools_FDR:
    input:
        expand(rules.qtltools_gather.output,_pass='permutations',allow_missing=True)
    output:
        'QTL/{QTL}/{tissue}_{variants}/permutations_all.{MAF}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    envmodules:
        'gcc/8.2.0',
        'r/4.2.2'
    shell:
        '''
        Rscript /cluster/work/pausch/alex/software/qtltools/scripts/qtltools_runFDR_cis.R {input} 0.05 {params.out}
        '''
