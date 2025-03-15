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
        TRF = config['reference_TRF'],
        RM = config['repeat_masking']
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
        expand(rules.happy.output['others'][2],sample=samples,chromosome=main_regions+['MT'],allow_missing=True)
    output:
        'happy/{read1}_{read2}.{filtering}.F1.csv'
    localrule: True
    shell:
        '''
        echo -e "variant size region truth query recall precision truth_TiTv query_TiTv F1_score sample chromosome" > {output}
        for i in {input}
        do
          awk -v I=$(basename $i) -F',' '$4=="PASS" {{ split(I,a,"."); print $1,$2,$3,$17,$38,$8,$9,$22,$43,$11,a[1],a[2] }}' $i >> {output}
        done
        '''

rule SV_comparison:
    input:
        panel = 'PanGenie_comparison/pangenome_panel.SVs.vcf',
        direct = 'PanGenie_comparison/HiFi_cohort.SVs.vcf'
    output:
        vcf = 'PanGenie_comparison/overlap.vcf',
        isec = 'PanGenie_comparison/overlap.isec'
    conda: 'jasmine'
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '60m'
    shell:
        '''
        jasmine --comma_filelist file_list={input.direct},{input.panel} threads={threads} out_file=/dev/stdout out_dir=$TMPDIR --keep_var_ids \
        genome_file={config[reference]} --pre_normalize --ignore_strand --ignore_type \
        max_dist_linear=0.5 max_dist=250 |\
        tee {output.vcf} |\
        grep -hoP "SUPP_VEC=\K\d+" |\
        awk ' {{ A[$1]+=1 }} END {{ print A["01"],A["10"],A["11"] }}' > {output.isec}
        '''

rule SV_QTL_comparison:
    input:
        jasmine = rules.SV_comparison.output['vcf'],
        QTL = '/cluster/work/pausch/HiFi_QTL/QTL/{QTL}/117_samples/testing_imputation/Testis_filtered/conditional.all.txt'
    output:
        pop_only_SVs = 'PanGenie_comparison/unique_SVs.{QTL}.list',
        significant_SVs = 'PanGenie_comparison/unique_SVs.{QTL}.significant.list',
    shell:
        '''
        awk '$8~/SUPP_VEC=10/ {{print $3}}' {input.jasmine} > {output.pop_only_SVs}
        LC_ALL=C; grep -Ff {output.pop_only_SVs} {input.QTL} | grep "1 1$" | sort -u > {output.significant_SVs}
        '''

rule find_disagreement_sites:
    input:
        expand('happy/{sample}.{chromosome}.mm2_bwa.bcf',sample=samples,chromosome=range(1,30))
    output:
        sites = 'Illumina_HiFi_comparison/disconcordant.sites'
    threads: 8
    resources:
        mem_mb = 2500,
        walltime = '24h'
    shell:
        '''
        bcf_printer() {{ bcftools view -m 2 -M 2 -v snps -e 'FORMAT/BD=="TP"' $1 | bcftools query -f '%CHROM %POS %INFO/Regions [%BVT ]' ; }}
        export -f bcf_printer

        parallel -j {threads} --line-buffer bcf_printer ::: {input} > {output}
        '''

rule group_SNP_sites:
    input:
        rules.find_disagreement_sites.output
    output:
        'Illumina_HiFi_comparison/disconcordant.bed'
    threads: 1
    resources:
        mem_mb = 2500
    shell:
        '''
        awk -v OFS='\\t' '{{print $1,$2,$2+1,$3,$4,$5}}' {input} | sort --parallel=2 -k1,1n -k2,2n | bedtools merge -d -1 -i /dev/stdin -o distinct,distinct,distinct,count -c 4,5,6,4 > {output}
        '''

rule SNP_disagreement:
    input:
        sites = rules.group_SNP_sites.output,
        QTL = '/cluster/work/pausch/HiFi_QTL/QTL/{QTL}/117_samples/testing_imputation/Testis_filtered/conditional.all.txt'
    output:
        HiFi_only_SNPs = '',
        significant_SNPs = ''
    shell:
        '''
        grep -v "," {input.sites} | awk '$7>60 {{print $1"_"$2"_SNP"}}' > {output.HiFi_only_SNPs}
        LC_ALL=C; grep -Ff {output.HiFi_only_SNPs} {input.QTL} | cut -d' ' -f 1 | sort | uniq -c | sort -k1,1nr > {output.significant_SNPs}
        '''''

rule temp:
    input:
        ''
    output:
        ''
    shell:
        '''
        awk '$16<1 {print $1}' /cluster/work/pausch/HiFi_QTL/QTL/eQTL/117_samples/testing_HiFi_illumina/Testis_HiFi_filtered/permutations_all.{1..29}.01.05.significant.txt | sort -u  > eGenes.HiFi.txt
        awk '$16<1e-10 {print $1}' /cluster/work/pausch/HiFi_QTL/QTL/eQTL/117_samples/testing_HiFi_illumina/Testis_HiFi_filtered/permutations_all.{1..29}.01.05.significant.txt | sort -u  > eGenes.HiFi.strict.txt

        awk '$16<1 {print $1}' /cluster/work/pausch/HiFi_QTL/QTL/eQTL/117_samples/testing_HiFi_illumina/Testis_Illumina_filtered/permutations_all.{1..29}.01.05.significant.txt | sort -u  > eGenes.Illumina.txt
        awk '$16<1e-10 {print $1}' /cluster/work/pausch/HiFi_QTL/QTL/eQTL/117_samples/testing_HiFi_illumina/Testis_Illumina_filtered/permutations_all.{1..29}.01.05.significant.txt | sort -u  > eGenes.Illumina.strict.txt
        '''
