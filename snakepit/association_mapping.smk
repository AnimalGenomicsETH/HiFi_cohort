from pathlib import PurePath

wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    tissue = r'Testis',
    variants = r'imputed|filtered',
    chromosome = r'\d+|X|Y',
    MAF = r'\d+',
    vcf = r'(eQTL|gwas)/\S+'

regions = list(map(str,range(1,30))) + ['X','Y']

rule all:
    input:
        expand('QTL/eQTL/Testis_{variants}_{covariates}/{_pass}.{chromosome}.{MAF}.{FDR}.txt.gz',_pass=('nominals',),chromosome=regions,variants=config['variants'],MAF=config['MAF'],FDR=('None',),covariates=config['covariates']['eQTL']['Testis']),
        expand('QTL/eQTL/Testis_{variants}_{covariates}/{_pass}.{chromosome}.{MAF}.{FDR}.txt.gz',_pass=('conditionals',),chromosome=regions,variants=config['variants'],MAF=config['MAF'],FDR=config['FDR'],covariates=config['covariates']['eQTL']['Testis'])

rule normalise_vcf:
    input:
        lambda wildcards: expand(config['variants'][wildcards.variants],**wildcards,allow_missing=True)
    output:
        multiext('QTL/variants/{variants}.{chromosome}.vcf.gz','','.tbi')
    threads: 4
    resources:
        mem_mb = 2500,
        disk_scratch = 50
    shell:
        '''
        bcftools norm --threads {threads} -f {config[reference]} -d exact -m -any {input} -Ou |\
        bcftools sort -T $TMPDIR -Ou - |\
        bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' -o {output[0]} -
        tabix -fp vcf {output[0]}
        '''

rule bcftools_set_id_for_plink:
    input:
        rules.normalise_vcf.output
    output:
         vcf = multiext('QTL/variants/{variants}.{chromosome}.hexID.vcf.gz','','.tbi'),
         ID_map = 'QTL/variants/{variants}.{chromosome}.ID_map.txt'
    threads: 2
    resources:
        mem_mb = 1500
    shell:
        '''
        bcftools annotate --threads {threads} --set-id '%VKX' -o {output.vcf[0]} {input[0]}
        paste <(bcftools query -f '%ID' {input[0]}) <(bcftools query -f '%ID' {output.vcf[0]}) > {output.ID_map}
        '''

rule exclude_MAF:
    input:
        rules.normalise_vcf.output
    output:
        'QTL/variants/{variants}.exclude_sites.{chromosome}.{MAF}.txt'
    shell:
        '''
        bcftools view --threads 2 -Q 0.{wildcards.MAF}:minor -Ou {input[0]} |\
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
        vcf = rules.normalise_vcf.output,
        exclude = rules.exclude_MAF.output,
        bed = lambda wildcards: expand(config['mol_QTLs'][wildcards.QTL][wildcards.tissue],**wildcards,allow_missing=True),
        cov = lambda wildcards: expand(config['covariates'][wildcards.QTL][wildcards.tissue][wildcards.covariates],**wildcards,allow_missing=True),
        mapping = lambda wildcards: 'QTL/{QTL}/{tissue}_{variants}_{covariates}/permutations_all.{chromosome}.{MAF}.{FDR}.thresholds.txt' if wildcards._pass == 'conditionals' else []
    output:
        merged = 'QTL/{QTL}/{tissue}_{variants}_{covariates}/{_pass}.{chromosome}.{MAF}.{FDR}.txt.gz'
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input),
        grp = lambda wildcards: '--grp-best' if wildcards.QTL == 'sQTL' else ''
    threads: 1
    resources:
        mem_mb = 2500,
    shell:
        '''
        QTLtools cis --vcf {input.vcf[0]} --exclude-sites {input.exclude} --bed {input.bed} --cov {input.cov} --std-err {params._pass} {params.grp} --window {config[window]} --normal --silent --log /dev/stderr --out /dev/stdout | pigz -p 2 -c > {output}
        '''

rule qtltools_FDR:
    input:
        expand(rules.qtltools_parallel.output,_pass='permutations',FDR='None',allow_missing=True)
    output:
        'QTL/{QTL}/{tissue}_{variants}_{covariates}/permutations_all.{chromosome}.{MAF}.{FDR}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    #conda: 'R'
    envmodules:
        'gcc/8.2.0',
        'R/4.2.2'
    shell:
        '''
        Rscript {config[FDR_script]} {input} 0.{wildcards.FDR} {params.out}
        '''
