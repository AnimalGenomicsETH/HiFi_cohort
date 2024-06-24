import pandas as pd
from pathlib import PurePath

HiFi_dataframe = pd.read_csv(config['sample_sheet'])
samples = sorted(list(HiFi_dataframe['Sample'].unique()))

regions = list(map(str,range(1,30))) + ['X','Y','MT','unplaced']
main_regions = list(map(str,range(1,30))) + ['X','Y_PAR','Y_HAPLOID']

include: 'snakepit/alignment.smk'
include: 'snakepit/variant_calling.smk'
include: 'snakepit/read_comparison.smk'
include: 'snakepit/methylation.smk'
include: 'snakepit/coverage.smk'
include: 'snakepit/eVariant_F1_overlap.smk'

#workflow._singularity_args = f'-B $TMPDIR -B {PurePath(config["reference"]).parent}'
workflow._singularity_args =  '-B $TMPDIR -B /cluster/work/pausch/inputs'

wildcard_constraints:
    sample = r'BSWCHEM\d+',
    regions = r'|'.join(regions) + r'|Y_PAR|Y_HAPLOID',
    mapper = r'bwa|strobe|mm2|pbmm2|wm2'

rule all:
    input:
        expand(rules.short_read_align.output,sample=samples,mapper=['bwa','strobe']),
        expand(rules.samtools_merge.output,sample=samples,mapper=['mm2','pbmm2'])
        #'variants/cohort.merged.vcf.gz',
        #'happy/F1.csv',
        #'methylation/merged.bed'
