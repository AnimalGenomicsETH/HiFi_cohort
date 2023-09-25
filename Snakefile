
import pandas as pd
from pathlib import PurePath

HiFi_dataframe = pd.read_csv(config['sample_sheet'])
samples = list(HiFi_dataframe['Sample'].unique())

include: 'snakepit/alignment.smk'
include: 'snakepit/variant_calling.smk'
include: 'snakepit/methylation.smk'

workflow._singularity_args = f'-B $TMPDIR -B {PurePath(config["reference"]).parent}'

rule all:
    input:
        expand(rules.samtools_merge.output,sample=samples),
        'variants/cohort.merged.vcf.gz',
        'happy/F1.csv',
        'methylation/merged.bed'
