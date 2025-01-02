ruleorder: bcftools_filter > sniffles_merge > sniffles_call
ruleorder: bcftools_concat > beagle4_impute

def get_DV_input(wildcards):
    match wildcards.mapper:
        case 'pbmm2' | 'mm2':
            return expand(rules.samtools_merge.output,sample=samples,mapper=wildcards.mapper)
        case 'bwa' | 'strobe':
            return expand(rules.short_read_align.output,sample=samples,mapper=wildcards.mapper)

rule deepvariant:
    input:
        get_DV_input
    output:
        expand('{mapper}_DV/{region}.Unrevised.vcf.gz',region=regions,allow_missing=True)
    params:
        name = lambda wildcards, output: PurePath(output[0]).parent,
        model = lambda wildcards: 'WGS' if wildcards.mapper in ['bwa','strobe'] else 'PACBIO',
        bam = lambda wildcards, input: PurePath(input[0]).suffix,
        index = lambda wildcards, input: PurePath(input[len(samples)]).suffix,
        config = 'config/deepvariant.yaml'
    localrule: True
    shell:
        '''
        snakemake -s /cluster/work/pausch/alex/BSW_analysis/snakepit/deepvariant.smk --configfile {params.config} \
        --config Run_name="{params.name}" model="{params.model}" \
        bam_name="{{sample}}.{wildcards.mapper}{params.bam}" bam_index="{params.index}" \
        --profile "slurm/fullNT" --resources storage_load=500 --nolock
        '''

rule bcftools_filter_DV:
    input:
        '{mapper}_DV/{region}.Unrevised.vcf.gz'
    output:
        '{mapper}_DV/{region}.filtered.vcf.gz'
    params:
        filter_expr = lambda wildcards: "'F_MISSING<0.2'"
    shell:
        '''
        bcftools view -i {params.filter_expr} -o {output} {input}
        tabix -p vcf {output}
        '''

rule bcftools_split_Y_PAR:
    input:
        '{mapper}_DV/Y.Unrevised.vcf.gz'
    output:
        expand('{{mapper}}_DV/{region}.Unrevised.vcf.gz',region=('Y_PAR','Y_HAPLOID'))
    shell:
        '''
        bcftools +scatter {input} -o $TMPDIR -Oz -S <(echo -e "Y:1-6822380\\tY_PAR\\nY:6822380-59476289\\tY_HAPLOID") --write-index
        mv $TMPDIR/Y_PAR.vcf.gz {output[0]}
        mv $TMPDIR/Y_HAPLOID.vcf.gz {output[1]}
        '''

rule beagle4_impute:
    input:
        '{mapper}_DV/{region}.Unrevised.vcf.gz'
    output:
        multiext('{mapper}_DV/{region}.beagle4.vcf.gz','','.tbi')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        name = lambda wildcards, output: PurePath(output[0]).name
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = '4h'
    shell:
        '''
        java -jar -Xss25m -Xmx50G /cluster/work/pausch/alex/software/beagle.27Jan18.7e1.jar ne=100 gl={input} nthreads={threads} out={params.prefix}
        cp {output[0]} $TMPDIR/{params.name}
        bcftools reheader -f {config[reference]}.fai -o {output[0]} $TMPDIR/{params.name}
        tabix -fp vcf {output[0]}
        '''

rule bcftools_filter_beagle:
    input:
        rules.beagle4_impute.output
    output:
        multiext('{mapper}_DV/{region}.beagle4_filtered.vcf.gz','','.tbi')
    params:
        filter_expr = lambda wildcards: f"-e 'DR2<0.7'{' -g ^het' if wildcards.region in ('X','Y_HAPLOID') else ''}"
    shell:
        '''
        bcftools view {params.filter_expr} -o {output[0]} {input[0]}
        tabix -p vcf {output[0]}
        '''

rule bcftools_concat:
    input:
        expand(rules.beagle4_impute.output[0],region=regions,allow_missing=True)
    output:
        multiext('{mapper}_DV/all.beagle4.vcf.gz','','.tbi')
    localrule: True
    shell:
        '''
        bcftools concat --threads {threads} --naive-force --no-version {input} > {output[0]}
        tabix -p vcf {output[0]}
        '''

rule pbsv_discover:
    input:
        expand(rules.samtools_merge.output,mapper='pbmm2',allow_missing=True)
    output:
        'SVs/{sample}.svsig.gz'
    conda:
        'pbccs'
    threads: 1
    resources:
        mem_mb = 15000
    shell:
        '''
        pbsv discover --ccs {input[0]} {output}
        '''

rule pbsv_call:
    input:
        signatures = expand(rules.pbsv_discover.output,sample=samples)
    output:
        'SVs/cohort.pbsv.vcf'
    conda:
        'pbccs'
    threads: 4
    resources:
        mem_mb = 40000,
        walltime = '4h'
    shell:
        '''
        pbsv call --hifi -j {threads} --log-level INFO --max-ins-length 20k {config[reference]} {input.signatures} {output}
        '''

rule sniffles_call:
    input:
        bam = expand(rules.samtools_merge.output,allow_missing=True),
        TR = 'GCA_002263795.4_ARS-UCD2.0_genomic.trf.bed'
    output:
        vcf = '{mapper}_SVs/{sample}.sniffles.denovo.vcf.gz',
        snf = '{mapper}_SVs/{sample}.sniffles.snf'
    threads: 4
    resources:
        mem_mb = 2000
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.bam[0]} --reference {config[reference]} --tandem-repeats {input.TR} --sample-id {wildcards.sample} --threads {threads} --max-del-seq-len 100000 --vcf {output.vcf} --snf {output.snf}
        '''

rule sniffles_merge:
    input:
        snfs = expand(rules.sniffles_call.output['snf'],sample=samples,allow_missing=True),
        TR = 'GCA_002263795.4_ARS-UCD2.0_genomic.trf.bed'
    output:
        vcf = '{mapper}_SVs/cohort.sniffles.denovo.vcf.gz'
    threads: 12
    resources:
        mem_mb = 3000
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.snfs} --reference {config[reference]} --tandem-repeats {input.TR} --threads {threads} --max-del-seq-len 100000 --vcf {output.vcf}
        '''

rule sniffles_filter:
    input:
        rules.sniffles_merge.output[0]
    output:
        '{mapper}_SVs/cohort.sniffles.denovo.filtered.vcf.gz'
    threads: 1
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools +fill-from-fasta {input} -- -c REF -f {config[reference]} |\
        bcftools view --threads {threads} -i 'abs(INFO/SVLEN)<=1000000&&INFO/SVTYPE!="BND"' -o {output} --write-index
        '''

rule sniffles_genotype:
    input:
        bam = expand(rules.samtools_merge.output,allow_missing=True),
        TR = config['reference_TRF'],
        SV_panel = rules.sniffles_filter.output
    output:
        vcf = '{mapper}_SVs/{sample}.sniffles.forced.vcf.gz'
    threads: 4
    resources:
        mem_mb = 2500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.bam[0]} --reference {config[reference]} --tandem-repeats {input.TR} --sample-id {wildcards.sample} --threads {threads} --max-del-seq-len 100000 --genotype-vcf {input.SV_panel} --vcf {output.vcf}
        '''

rule bcftools_merge_sniffles:
    input:
        expand(rules.sniffles_genotype.output['vcf'],sample=samples,allow_missing=True)
    output:
        vcf = '{mapper}_SVs/cohort.sniffles.forced.vcf.gz'
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools merge --write-index -o {output} {input}
        '''

rule bcftools_filter:
    input:
        rules.bcftools_merge_sniffles.output
    output:
        multiext('{mapper}_SVs/filtered/{region}.vcf.gz','','.csi')
    params:
        regions = ','.join(regions[:-1]),
        _dir = lambda wildcards, output: PurePath(output[0]).parent#ith_suffix('').with_suffix('').with_suffix('').with_suffix('')
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools +fill-from-fasta {input} -- -c REF -f {config[reference]} |\
        bcftools view --threads {threads} -i 'abs(INFO/SVLEN)<=100000&&INFO/SVTYPE!="BND"' |\
        bcftools +scatter - -o {params._dir} -Oz --threads {threads} --write-index -s {params.regions} -x unplaced --no-version
        '''

rule merge_QTL_variants:
    input:
        small = rules.beagle4_impute.output,
        SV = rules.bcftools_filter.output
    output:
        multiext('QTL_variants/{region}.{mapper}.merged.vcf.gz','','.tbi')
    threads: 2
    resources:
        mem_mb = 1500,
        walltime = '30m'
    shell:
        '''
        bcftools norm --threads {threads} -f {config[reference]} -m -any -Ou {input.small[0]} |\
        bcftools sort -T $TMPDIR -Ou - |\
        bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' -o $TMPDIR/small.vcf.gz
        tabix -p vcf $TMPDIR/small.vcf.gz

        bcftools concat --allow-overlaps --threads {threads} -Ou {input.SV[0]} $TMPDIR/small.vcf.gz |\
        bcftools +fill-tags - -o {output[0]} -- -t all

        tabix -p vcf {output[0]}
        '''

rule bcftools_concat_QTL:
    input:
        expand(rules.merge_QTL_variants.output[0],mapper='mm2',region=regions[:-1])
    output:
        multiext('QTL_variants/chromosomes.vcf.gz','','.tbi')
    localrule: True
    shell:
        '''
        bcftools concat --threads {threads} --naive-force --no-version {input} > {output[0]}
        tabix -p vcf {output[0]}
        '''

rule split_SV_sequences:
    input:
        expand(rules.bcftools_merge_sniffles.output,mapper='mm2',allow_missing=True)
    output:
        'mm2_SVs/SV_sequences/{region}.fa'
    shell:
        '''
        bcftools query -f '%ID\\t%REF\\t%ALT\\t%INFO/SVTYPE\\t%INFO/SVLEN' {input} |\
        awk '$4!="."&&$3!~/</ {{if (length($2)>length($3)) {{print ">"$1"_"$4"_"$5"\\n"$2}} else {{print ">"$1"_"$4"_"$5"\n"$3}} }}' |\
        seqtk seq -U - > {output}
        '''

rule repeat_masker:
    input:
        rules.split_SV_sequences.output
    output:
        '{mapper}_SVs/SV_sequences/{region}.fa.out'
    threads: 2
    resources:
        mem_mb = 1500,
        walltime = '4h'
    shell:
        '''
        RepeatMasker -pa $(({threads}/1)) -e rmblast -lib {config[repeat_library]} -q -no_is {input}
        if [ ! -f {output} ]; then
          seqtk seq -l60 {input} > {output}
        fi
        '''

rule TRF:
    input:
        rules.split_SV_sequences.output
    output:
        '{mapper}_SVs/SV_sequences/{region}.fa.trf'
    shell:
        '''
        TRF {input} 2 5 7 80 10 50 2000 -h -ngs > {output}
        '''

rule merge_masked_chromosomes:
    input:
        out = expand(rules.repeat_masker.output,region=regions[:-2],allow_missing=True),
        trf = expand(rules.TRF.output,region=regions[:-2],allow_missing=True),
    output:
        repeats = '{mapper}_SVs/SV_sequences.repeats.gz',
        VNTRs = '{mapper}_SVs/SV_sequences.VNTRs.gz'
    localrule: True
    shell:
        '''
        echo "SV length element class" > {output.repeats}
        
        awk 'NR>3 {{print $5,$7-$6,$10,$11}}' {input.out} >> {output.repeats}

        echo "SV period L score motif" > {output.VNTRs}

        awk '{{if ($1~/@/) {{S=substr($1,2); next}} {{print S,$3,$4,$8,$14}} }}' {input.trf} >> {output.VNTRs}
        '''

rule split_chromosomes:
    input:
        rules.samtools_merge.output
    output:
        expand('phasing/splits/{{sample}}.{{mapper,mm2|pbmm2}}.{region}.{ext}',ext=('bam','bam.csi'),region=regions[:-1])
    params:
        regions = regions[:-1]
    threads: 2
    resources:
        mem_mb = 2500,
        walltime = '2h'
    shell:
        '''
        for i in {params.regions}
        do
          samtools view -@ {threads} --write-index -o phasing/splits/{wildcards.sample}.{wildcards.mapper}.${{i}}.bam {input} $i
        done
        '''

rule hiphase_chromosome:
    input:
        #bam = expand(rules.samtools_merge.output,mapper='mm2',allow_missing=True),
        bam = expand('phasing/splits/{sample}.mm2.{{region}}.{ext}',ext=('bam','bam.csi'),sample=samples),#expand(rules.split_chromosomes.output,mapper='mm2',sample=samples,allow_missing=True),
        vcf = expand(rules.merge_QTL_variants.output,mapper='mm2',allow_missing=True)
    output:
        bam = expand('phasing/{sample}.{{region}}.mm2.phased.{ext}',ext=('bam','bam.bai'),sample=samples),
        vcf = expand("phasing/cohort.{{region}}.mm2.phased.{ext}",ext=('vcf.gz','vcf.gz.tbi')),
        #summary = multiext("phasing/{sample}.{region}.mm2.phased",".summary.tsv",".stats.tsv",".blocks.tsv")
    params:
        IO_list = lambda wildcards, input, output: ' '.join(f'--bam {input.bam[I]} --output-bam {output.bam[I]} --sample-name {S}' for I,S in enumerate(samples))
    threads: 12
    resources:
        mem_mb = 5000,
        walltime = '24h'
    shell:
        '''
        hiphase --threads {threads} --io-threads 6 --reference {config[reference]} \
        {params.IO_list} \
        --vcf {input.vcf[0]} --output-vcf {output.vcf[0]} \
        --global-realignment-cputime 300 \
        --min-vcf-qual 20 \
        --phase-singletons \
        --summary-file phasing/summary.{wildcards.region}.txt --stats-file phasing/stats.{wildcards.region}.txt --blocks-file phasing/blocks.{wildcards.region}.txt
        '''

rule summarise_phasing_chromosome:
    input:
        expand('phasing/summary.{region}.txt',region=regions[:-1])
    output:
        'phasing/summary.csv'
    localrule: True
    shell:
        '''
        echo "sample_name chromosome num_variants num_heterozygous num_phased num_blocks basepairs_per_block_max block_ng50" > {output}
        for R in {input}
        do
          filename=$(basename -- "$R")
          filename="${{filename%.*}}"
          filename="${{filename##*.}}"
          awk -v I=$filename '$2==I {{print $1,$2,$3,$4,$5,$9,$19,$NF}}' $R  >> {output}
        done
        '''

rule concat_phases:
    input:
        expand(rules.hiphase_chromosome.output.vcf[0],region=regions[:-1])
    output:
        'phasing/cohort.chromosomes.phased.vcf.gz'
    localrule: True
    shell:
        '''
        bcftools concat --naive-force -o {output} {input}
        '''

rule samtools_merge_phased:
    input:
        expand('phasing/{{sample}}.{region}.mm2.phased.bam',region=regions[:-1])
    output:
        multiext('phasing/{sample}.mm2.phased.bam','','.csi')
    threads: 4
    resources:
        mem_mb = 2500
    shell:
        '''
        samtools merge -@ {threads} --write-index -o {output[0]} {input}
        '''

rule hiphase:
    input:
        bam = expand(rules.samtools_merge.output,mapper='mm2',allow_missing=True),
        small = rules.bcftools_concat_QTL.output, 
    output:
        bam = multiext("alignments/{sample}.mm2.haplotagged.bam",'','.bai'),
        small = multiext("variants/{sample}.mm2.DV.phased.vcf.gz",'','.tbi'),
        summary = multiext("variants/{sample}.mm2.phased",".summary.tsv",".stats.tsv",".blocks.tsv")
    threads: 8
    resources:
        mem_mb = 4000,
        walltime = '8h'
    shell:
        '''
        hiphase --threads {threads} --io-threads 6 --reference {config[reference]} \
        --bam {input.bam[0]} --output-bam {output.bam[0]} \
        --vcf {input.small[0]} --output-vcf {output.small[0]} \
        --sample-name {wildcards.sample} \
        --global-realignment-cputime 300 \
        --min-vcf-qual 20 \
        --phase-singletons \
        --summary-file {output.summary[0]} --stats-file {output.summary[1]} --blocks-file {output.summary[2]}
        '''
