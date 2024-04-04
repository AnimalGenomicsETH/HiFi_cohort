from pathlib import PurePath

rule bcftools_split:
    input:
        '{mapper}_DV/{chromosome}.{filtering}.vcf.gz'
    output:
        expand('{mapper}_DV/PER_SAMPLE_{chromosome}_{filtering}/{sample}.vcf.gz',sample=samples,allow_missing=True)
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools +split -i 'GT[*]="alt"' -Oz -o {params._dir} {input}
        '''

rule make_happy_regions:
    input:
        TRF = '/cluster/work/pausch/alex/Pop_HiFi/GCA_002263795.4_ARS-UCD2.0_genomic.trf.bed',
        RM = '/cluster/work/pausch/alex/Pop_HiFi/GCF_002263795.3.repeatMasker.out.gz'
    output:
        regions = 'happy/regions.tsv',
        beds = expand('happy/{region}.bed',region=('VNTR','TE','normal'))
    localrule: True
    shell:
        '''
        bedtools sort -faidx {config[reference]}.fai -i {input.TRF} > {output.beds[0]}

        zgrep -P "(LTR|LINE|SINE)" {input.RM} |\
        RM2Bed.py - - |\
        bedtools sort -faidx {config[reference]}.fai -i - |\
        bedtools subtract -a /dev/stdin -b {output.beds[0]} -A |\
        cut -f -3 > {output.beds[1]}

        cat {output.beds[0]} {output.beds[1]} |\
        bedtools sort -faidx {config[reference]}.fai -i - |\
        bedtools complement -g {config[reference]}.fai -i /dev/stdin > {output.beds[2]}

        echo -e "VNTR\\t{output.beds[0]}\\nTE\\t{output.beds[1]}\\nnormal\\t{output.beds[2]}" > {output.regions}

        '''

rule happy:
    input:
        vcf1 = '{read1}_DV/PER_SAMPLE_{chromosome}_{filtering}/{sample}.vcf.gz',
        vcf2 = '{read2}_DV/PER_SAMPLE_{chromosome}_{filtering}/{sample}.vcf.gz',
        reference = config['reference_uncompressed'],
        regions = 'happy/regions.tsv'
    output:
        csv = 'happy/{sample}.{chromosome}.{filtering}.{read1}_{read2}.summary.csv',
        others = multiext('happy/{sample}.{chromosome}.{filtering}.{read1}_{read2}','.bcf','.bcf.csi','.extended.csv','.roc.all.csv.gz','.runinfo.json')
    params:
        _dir = lambda wildcards, output: PurePath(output.csv).with_suffix('').with_suffix('')
    container: '/cluster/work/pausch/alex/software/images/hap.py_latest.sif'
    threads: 2
    resources:
        mem_mb = 15000,
        scratch = '10G',
        walltime = '30m',
        storage_load = 1
    shell:
        '''
        /opt/hap.py/bin/hap.py -r {input.reference} --bcf --usefiltered-truth --no-roc --no-json -L --pass-only --scratch-prefix $TMPDIR -X --threads {threads} --stratification {input.regions} -o {params._dir} {input.vcf1} {input.vcf2}
        '''

rule gather_happy:
    input:
        expand(rules.happy.output['others'][2],sample=samples,chromosome=main_regions,allow_missing=True)
    output:
        'happy/{read1}_{read2}.{filtering}.F1.csv'
    localrule: True
    shell:
        '''
        echo -e "variant region truth query recall precision truth_TiTv query_TiTv F1_score sample chromosome" > {output}
        for i in {input}
        do
          awk -v I=$(basename $i) -F',' '$2=="*"&&$4=="PASS" {{ split(I,a,"."); print $1,$3,$17,$38,$8,$9,$22,$43,$11,a[1],a[2] }}' $i >> {output}
        done
        '''
