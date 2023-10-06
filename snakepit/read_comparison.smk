from pathlib import PurePath


rule bcftools_split:
    input:
        '{mapper}_DV/{chromsome}.{imputed}.vcf.gz'
    output:
        expand('{mapper}_DV/PER_SAMPLE_{chromosome}/{sample}.vcf.gz',sample=config['samples'],allow_missing=True)
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
        scratch = '10G'
    shell:
        '''
        /opt/hap.py/bin/hap.py -r {input.reference} --bcf --usefiltered-truth --no-roc --no-json -L --pass-only --scratch-prefix $TMPDIR -X --threads {threads} -o {params._dir} {input.vcf1} {input.vcf2}
        '''

rule gather_happy:
    input:
        expand(rules.happy.output[0],sample=samples,chromosomes=regions)
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

rule bcftools_isec:
    input:
        'A'
    output:
        'B'
    shell:
        '''
        bcftools isec -n +1 -c some {input.HiFi} {input.SR} > {output[0]}
        '''

rule samtools_coverage:
    input:
        lambda wildcards: expand('/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long-read_alignments/PacBio_CCS/eQTL/{sample}.mm2.cram' if wildcards.read == 'HiFi' else '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD_eQTL/{sample}.bam',sample=config['samples'])
    output:
        'coverage/{read}.txt'
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '24h'
    shell:
        '''
        samtools coverage --reference {config[reference]} -q 20 -b <(echo {input} | tr ' ' '\\n') > {output}
        '''

rule prep_cov:
    output:
        'coverage/regions.bed'
    localrule: True
    shell:
        '''
        awk '/^[1-9]/||/^X/||/^Y/ {{ print $1"\t0\t"$2 }} ' {config[reference]}.fai > {output}
        '''

rule samtools_bedcov:
    input:
        bed = rules.prep_cov.output,
        bam = lambda wildcards: '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long-read_alignments/PacBio_CCS/eQTL/{sample}.mm2.cram' if wildcards.read == 'HiFi' else '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD_eQTL/{sample}.bam'
    output:
        'coverage/{sample}.{read}.{size}.bedcov'
    threads: 1
    resources:
        mem_mb = 2500,
        walltime = '30m'
    shell:
        '''
        samtools bedcov --reference {config[reference]} -Q 20 -d {wildcards.size} {input.bed} {input.bam} | awk ' {{ print "{wildcards.sample}",$1,$5/$3 }} ' > {output}
        '''

rule gather_bedcovs:
    input:
        hifi = expand(rules.samtools_bedcov.output,sample=config['samples'],read='HiFi',allow_missing=True),
        SR = expand(rules.samtools_bedcov.output,sample=config['samples'],read='SR',allow_missing=True)
    output:
        'coverage/total.{size}.tsv'
    localrule:  True
    shell:
        '''
        echo "read size sample chromosome covered" > {output}
        awk ' {{ print "HiFi","{wildcards.size}",$0 }} ' {input.hifi} >> {output}
        awk ' {{ print "SR","{wildcards.size}",$0 }} ' {input.SR} >> {output}
        '''
