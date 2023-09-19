
include: 'snakepit/methylation.smk'
include: 'snakepit/alignment.smk'
include: 'snakepit/variant_calling.smk'

HiFi_dataframe = pd.read_csv(config['sample_sheet'])
samples = list(HiFi_dataframe['samples'].unique())

workflow._singularity_args = f'-B $TMPDIR -B {PurePath(config["reference"]).parent}'

rule all:
    input:
        expand(rules.samtools_merge.output,samples=samples),
        'variants/cohort.merged.vcf.gz',
        'happy/F1.csv',
        'methylation/merged.bed'
