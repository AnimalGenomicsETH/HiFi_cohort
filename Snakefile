import pandas as pd
from pathlib import PurePath

HiFi_dataframe = pd.read_csv(config['sample_sheet'])
samples = sorted(list(HiFi_dataframe['Sample'].unique()))

regions = list(map(str,range(1,30))) + ['X','Y','MT','unplaced']
main_regions = list(map(str,range(1,30))) + ['X','Y_PAR','Y_HAPLOID']

wildcard_constraints:
    sample = r'BSWCHEM\d+',
    regions = r'|'.join(regions) + r'|Y_PAR|Y_HAPLOID',
    mapper = r'bwa|strobe|mm2|pbmm2|wm2',
    _pass = r'permutations|conditionals|nominals',
    tissue = r'Testis',
    variants = r'imputed|filtered',
    chromosome = r'\d+|X|Y\w*',
    MAF = r'\d+',
    vcf = r'(eQTL|gwas)/\S+'

include: 'snakepit/alignment.smk'
include: 'snakepit/variant_calling.smk'
include: 'snakepit/read_comparison.smk'
include: 'snakepit/coverage.smk'
include: 'snakepit/eVariant_F1_overlap.smk'
include: 'snakepit/VEP.smk'
include: 'snakepit/association_mapping.smk'

workflow._singularity_args =  '-B $TMPDIR'

rule all:
    input:
        ## Illumina and HiFi alignments
        expand(rules.short_read_align.output,sample=samples,mapper='bwa'),
        expand(rules.samtools_merge.output,sample=samples,mapper='mm2'),
        ## F1 comparison
        'happy/mm2_bwa.Unrevised.F1.csv',
        ## SV analysis
        'mm2_SVs/SV_sequences.VNTRs.gz',
        'VEP/genomic.INS.bed',
        ## eQTL and sQTL mapping
        expand( 'QTL/{QTL}/Testis_filtered_PCA/conditionals.{chromosome}.01.05.txt.gz',QTL=('eQTL','sQTL'),chromosome=regions[:-2])

