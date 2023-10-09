from pathlib import PurePath


rule bcftools_split:
    input:
        '{mapper}_DV/{chromosome}.beagle4.vcf.gz'
    output:
        expand('{mapper}_DV/PER_SAMPLE_{chromosome}/{sample}.vcf.gz',sample=samples,allow_missing=True)
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools +split -i 'GT[*]="alt"' -Oz -o {params._dir} {input}
        '''

rule happy:
    input:
        vcf1 = '{read1}_DV/PER_SAMPLE_{chromosome}/{sample}.vcf.gz',
        vcf2 = '{read2}_DV/PER_SAMPLE_{chromosome}/{sample}.vcf.gz',
        reference = config['reference']
    output:
        csv = 'happy_{read1}_{read2}/{sample}.{chromosome}.summary.csv',
        others = multiext('happy_{read1}_{read2}/{sample}.{chromosome}','.bcf','.bcf.csi','.extended.csv','.roc.all.csv.gz','.runinfo.json')
    params:
        _dir = lambda wildcards, output: PurePath(output.csv).with_suffix('').with_suffix('')
    container: '/cluster/work/pausch/alex/software/images/hap.py_latest.sif'
    threads: 2
    resources:
        mem_mb = 5000,
        scratch = '10G',
        walltime = '30m',
        storage_load = 1
    shell:
        '''
        /opt/hap.py/bin/hap.py -r {input.reference} --bcf --usefiltered-truth --no-roc --no-json -L --pass-only --scratch-prefix $TMPDIR -X --threads {threads} -o {params._dir} {input.vcf1} {input.vcf2}
        '''

rule gather_happy:
    input:
        expand(rules.happy.output[0],sample=samples,chromosome=regions[:-1],allow_missing=True)
    output:
        'happy/{read1}_{read2}.F1.csv'
    localrule: True
    shell:
        '''
        echo -e "variant truth query recall precision truth_TiTv query_TiTv sample" > {output}
        for i in {input}
        do
          awk -v I=$(basename $i) -F',' '$2=="PASS" {{ split(I,a,"."); print $1,$3,$6,$10,$11,$14,$15,a[1] }}' $i >> {output}
        done
        '''
