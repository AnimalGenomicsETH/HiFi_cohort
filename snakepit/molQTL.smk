from pathlib import PurePath

wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    tissue = r'Testis',
    chromosome = r'\d+|X\w*|Y\w*',
    MAF = r'\d+',
    vcf = r'(eQTL|gwas)/\S+'

regions = list(map(str,range(1,30))) + ['X_het_missing','Y_PAR','Y_HAP_het_missing']

rule all:
    input:
        expand('{QTL}/117_samples/testing_imputation/Testis_{variants}/{_pass}.{chromosome}.{MAF}.{FDR}.txt.gz',_pass=('nominals',),chromosome=regions,variants=config['variants'],QTL=config['mol_QTLs'],MAF=config['MAF'],FDR=('None',)),
        expand('{QTL}/117_samples/testing_imputation/Testis_{variants}/{_pass}.{chromosome}.{MAF}.{FDR}.txt.gz',_pass=('conditionals',),chromosome=regions,variants=config['variants'],QTL=config['mol_QTLs'],MAF=config['MAF'],FDR=config['FDR'])

rule normalise_vcf:
    input:
        lambda wildcards: expand(config['variants'][wildcards.variants],**wildcards,allow_missing=True)
    output:
        multiext('variants/117_samples/HiFi_SVs/{variants}.{chromosome}.vcf.gz','','.tbi')
    threads: 4
    envmodules:
    	'stack/2024-06',
    	#'gcc/8.2.0',
    	'htslib'
    resources:
        mem_mb = 2500,
        disk_scratch = 50
    shell:
        '''
        /path/to/bcftools norm --threads {threads} -f {config[reference]} -d exact -m -any {input} -Ou |\
        /path/to/bcftools sort -T $TMPDIR -Ou -|\
  	    /path/to/bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' -o {output[0]} -
        tabix -fp vcf {output[0]}
        '''

rule exclude_MAF:
    input:
        rules.normalise_vcf.output
    output:
        'variants/117_samples/HiFi_SVs/{variants}.exclude_sites.{chromosome}.{MAF}.txt'
    shell:
        '''
        /path/to/bcftools view -g ^miss -Ou {input[0]} |\
        /path/to/bcftools view --threads 2 -Q 0.{wildcards.MAF}:minor -Ou - |\
        /path/to/bcftools query -f '%ID\n'  - > {output}
        '''

def get_pass(_pass,input):
    if _pass == 'permutations':
        return f'--permute {config["permutations"]}'
    elif _pass == 'conditionals':
        return f'--mapping {input.mapping}'
    elif _pass == 'nominals':
        return f'--nominal {config.get("nominal_threshold")}'

rule qtltools_parallel:
    input:
        vcf = rules.normalise_vcf.output,
        exclude = rules.exclude_MAF.output,
        bed = lambda wildcards: expand(config['mol_QTLs'][wildcards.QTL][wildcards.tissue],**wildcards,allow_missing=True),
        cov = lambda wildcards: expand(config['covariates'][wildcards.QTL][wildcards.tissue],**wildcards,allow_missing=True),
        mapping = lambda wildcards: '{QTL}/117_samples/testing_imputation/{tissue}_{variants}/permutations_all.{chromosome}.{MAF}.{FDR}.thresholds.txt' if wildcards._pass == 'conditionals' else []
    output:
        merged = '{QTL}/117_samples/testing_imputation/{tissue}_{variants}/{_pass}.{chromosome}.{MAF}.{FDR}.txt.gz'
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input),
        grp = lambda wildcards: '--grp-best' if wildcards.QTL == 'sQTL' else ''
    threads: 1
    resources:
        mem_mb = 2000,
        walltime = "24h"
    shell:
        '''
        /path/to/qtltools/bin/QTLtools cis --vcf {input.vcf[0]} --exclude-sites {input.exclude} --bed {input.bed} --cov {input.cov} --std-err {params._pass} {params.grp} --window {config[window]} --silent --log /dev/stderr --out /dev/stdout | pigz -p 2 -c > {output}
        '''

rule qtltools_FDR:
    input:
        expand(rules.qtltools_parallel.output,_pass='permutations',FDR='None',allow_missing=True)
    output:
        '{QTL}/117_samples/testing_imputation/{tissue}_{variants}/permutations_all.{chromosome}.{MAF}.{FDR}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    envmodules:
        'stack/2024-06',
        'r/4.3.2'
    shell:
        '''
        Rscript {config[FDR_script]} {input} 0.{wildcards.FDR} {params.out}
        '''