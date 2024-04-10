from pathlib import PurePath

wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    tissue = r'Testis',
    chromosome = r'\d+|X|Y',
    MAF = r'\d+',
    vcf = r'(eQTL|gwas)/\S+'

regions = list(map(str,range(1,30))) + ['X','Y']

rule all:
    input:
        expand('QTL/eQTL/Testis_{variants}/{_pass}.{chromosome}.{MAF}.txt.gz',_pass=('nominals','conditionals'),chromosome=regions,variants=config['variants'],MAF=config['MAF'])

rule normalise_vcf:
    input:
        lambda wildcards: expand(config['variants'][wildcards.variants],**wildcards,allow_missing=True)
    output:
        'QTL/variants/{variants}.{chromosome}.vcf.gz'
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
        lambda wildcards: expand(config['variants'][wildcards.variants],**wildcards,allow_missing=True)
    output:
        'QTL/variants/{variants}.exclude_sites.{chromosome}.{MAF}.txt'
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
        return f'--nominal {config.get("nominal_threshold",0.05)}'

rule qtltools_parallel:
    input:
        vcf = lambda wildcards: expand(config['variants'][wildcards.variants],**wildcards,allow_missing=True),
        exclude = rules.exclude_MAF.output,
        bed = lambda wildcards: expand(config['mol_QTLs'][wildcards.QTL][wildcards.tissue],**wildcards,allow_missing=True),
        cov = lambda wildcards: expand(config['covariates'][wildcards.QTL][wildcards.tissue],**wildcards,allow_missing=True),
        mapping = lambda wildcards: 'QTL/{QTL}/{tissue}_{variants}/permutations_all.{chromosome}.{MAF}.thresholds.txt' if wildcards._pass == 'conditionals' else []
    output:
        merged = 'QTL/{QTL}/{tissue}_{variants}/{_pass}.{chromosome}.{MAF}.txt.gz'
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input),
        grp = lambda wildcards: '--grp-best' if wildcards.QTL == 'sQTL' else ''
    threads: 1
    resources:
        mem_mb = 2500,
    shell:
        '''
        QTLtools cis --vcf {input.vcf} --exclude-sites {input.exclude} --bed {input.bed} --cov {input.cov} --std-err {params._pass} {params.grp} --window {config[window]} --normal --silent --log /dev/stderr --out /dev/stdout | pigz -p 2 -c > {output}
        '''

rule qtltools_FDR:
    input:
        expand(rules.qtltools_parallel.output,_pass='permutations',allow_missing=True)
    output:
        'QTL/{QTL}/{tissue}_{variants}/permutations_all.{chromosome}.{MAF}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    conda: 'R'
    shell:
        '''
        Rscript /cluster/work/pausch/alex/software/qtltools/scripts/qtltools_runFDR_cis.R {input} 0.05 {params.out}
        '''
