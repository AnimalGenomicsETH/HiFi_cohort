ruleorder: bcftools_filter > sniffles_merge > sniffles_call

regions = list(map(str,range(1,30))) + ['X','Y','MT','unplaced']

def get_DV_input(wildcards):
    match wildcards.mapper:
        case 'pbmm2' | 'mm2':
            return expand(rules.samtools_merge.output,sample=samples,mapper=wildcards.mapper)
        case 'bwa' | 'strobe':
            return expand(rules.short_read_align.output,sample=samples,mapper=wildcards.mapper)

rule deepvariant:
    input:
        get_DV_input,
        #lambda wildcards: expand(rules.samtools_merge.output,sample=samples,mapper='mm2' if wildcards.mapper=='HiFi' else 'pbmm2'),
        config = 'config/DV.yaml'
    output:
        expand('{mapper}_DV/{region}.Unrevised.vcf.gz',region=regions,allow_missing=True)
    params:
        name = lambda wildcards, output: PurePath(output[0]).parent,
        model = lambda wildcards: 'WGS' if wildcards.mapper in ['bwa','strobe'] else 'PACBIO',
        bam = lambda wildcards, input: PurePath(input[0]).suffix,
        index = lambda wildcards, input: PurePath(input[len(samples)]).suffix #'.csi' if wildcards.mapper in ['mm2','pbmm2'] else 'crai'
    localrule: True
    shell:
        '''
        snakemake -s /cluster/work/pausch/alex/BSW_analysis/snakepit/deepvariant.smk --configfile {input.config} \
        --config Run_name="{params.name}" model="{params.model}" \
        bam_name="{{sample}}.{wildcards.mapper}{params.bam}" bam_index="{params.index}" \
        --profile "slurm/fullNT" --resources storage_load=500 --nolock
        '''

rule beagle4_impute:
    input:
        '{read_type}_DV/{region}.Unrevised.vcf.gz'
    output:
        multiext('{read_type}_DV/{region}.beagle4.vcf.gz','','.tbi')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        name = lambda wildcards, output: PurePath(output[0]).name
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = '4h'
    shell:
        '''
        java -jar -Xss25m -Xmx50G /cluster/work/pausch/alex/software/beagle.27Jan18.7e1.jar gl={input} nthreads={threads} out={params.prefix}
        cp {output[0]} $TMPDIR/{params.name}
        bcftools reheader -f {config[reference]}.fai -o {output[0]} $TMPDIR/{params.name}
        tabix -fp vcf {output[0]}
        '''

rule bcftools_concat:
    input:
        expand(rules.beagle4_impute.output[0],allow_missing=True)
    output:
        multiext('{read_type}_DV/all.beagle4.vcf.gz','','.tbi')
    shell:
        '''
        bcftools concat {input} > {output}
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
        TR = 'GCA_002263795.4_ARS-UCD1.3_genomic.trf.bed'
    output:
        vcf = '{mapper}_SVs/{sample}.sniffles.vcf.gz',
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
        TR = 'GCA_002263795.4_ARS-UCD1.3_genomic.trf.bed'
    output:
        vcf = '{mapper}_SVs/cohort.sniffles.vcf.gz'
    threads: 12
    resources:
        mem_mb = 3000
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.snfs} --reference {config[reference]} --tandem-repeats {input.TR} --threads {threads} --max-del-seq-len 100000 --vcf {output.vcf}
        '''

rule bcftools_filter:
    input:
        rules.sniffles_merge.output
    output:
        multiext('{mapper}_SVs/InDels.sniffles.vcf.gz','','.tbi')
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools +fill-from-fasta {input} -- -c REF -f {config[reference]} |\
        bcftools view --threads {threads} -i 'abs(INFO/SVLEN)<1000000&&INFO/SVTYPE!="BND"' -o {output[0]}
        tabix -p vcf {output[0]}
        '''

checkpoint split_SV_sequences:
    input:
        rules.bcftools_filter.output
    output:
        directory('{mapper}_SVs/SV_sequences')
    shell:
        '''
        mkdir -p {output}
        bcftools query -f '%ID\\t%REF\\t%ALT\\n' {input[0]} | awk 'length($2)>length($3) {{print ">"$1"\\n"$2;next}} {{print ">"$1"\\n"$3}}' |\
        split -l 10000 -d -a 4 --additional-suffix ".fa" - {output}/
        '''

rule repeat_masker:
    input:
        '{mapper}_SVs/SV_sequences/{chunk}.fa',
    output:
        '{mapper}_SVs/SV_sequences/{chunk}.fa.out'
    threads: 2
    resources:
        mem_mb = 1500,
        walltime = '4h'
    shell:
        '''
        RepeatMasker -xsmall -pa $(({threads}/1)) -e rmblast -lib {config[repeat_library]} -qq -no_is {input}
        if [ ! -f {output} ]; then
          seqtk seq -l60 {input} > {output}
        fi
        '''

rule TRF:
    input:
        '{mapper}_SVs/SV_sequences/{chunk}.fa'
    output:
        '{mapper}_SVs/SV_sequences/{chunk}.fa.trf'
    shell:
        '''
        TRF {input} 2 5 7 80 10 50 2000 -h -ngs > {output}
        '''

def aggregate_out_chunks(wildcards):
    checkpoint_output = checkpoints.split_SV_sequences.get(**wildcards).output[0]
    return expand('{fpath}/{chunk}.fa.out',fpath=checkpoint_output,chunk=glob_wildcards(PurePath(checkpoint_output).joinpath('{chunk}.fa')).chunk)

def aggregate_trf_chunks(wildcards):
    checkpoint_output = checkpoints.split_SV_sequences.get(**wildcards).output[0]
    return expand('{fpath}/{chunk}.fa.trf',fpath=checkpoint_output,chunk=glob_wildcards(PurePath(checkpoint_output).joinpath('{chunk}.fa')).chunk)

rule merge_masked_chromosomes:
    input:
        out = aggregate_out_chunks,
        trf = aggregate_trf_chunks
    output:
        out = '{mapper}_SVs/SV_sequences.out.gz',
        trf = '{mapper}_SVs/SV_sequences.trf.gz'
    localrule: True
    shell:
        '''
        cat {input.out} | bgzip --threads {threads} -c > {output.out}
        cat {input.trf} | bgzip --threads {threads} -c > {output.trf}

        #awk -v c=0 '{{if (/@/) {{S="N"; next}} }} {{ if (S=="N"&&$5>6&&$4>=5) {{++c; S="Y"}} }} END {{print c}}' *.fa.trf 
        #grep -hP "\sLTR/" *out | /cluster/work/pausch/alex/software/RepeatMasker/util/buildSummary.pl - | awk '/^Sniffles2/&&$3>100 {{print $1}}' > ltr
        #grep -hP "/RTE" *out | /cluster/work/pausch/alex/software/RepeatMasker/util/buildSummary.pl - | awk '/^Sniffles2/&&$3>100 {{print $1}}' > rte
        #comm -12 <(sort ltr) <(sort rte) | wc -l
        '''

rule hiphase:
    input:
        bam = rules.samtools_merge.output,
        vcf_snv = 'vcf.gz',
        reference = config['reference']
    output:
        cram = "{sample}/{sample}.mm2.haplotagged.cram",
        vcf_snv = "{sample}/{sample}.mm2.DV.phased.vcf.gz",
        summary = multiext("{sample}/{sample}.mm2.phased",".summary.tsv","stats.tvt","blocks.tsv")
    threads: 16
    resources:
        mem_mb = 1500,
        walltime = '4h'
    shell:
        '''
        hiphase -t {threads} -r {input.reference} -b {input.cram} -p {output.cram} -c {input.vcf} -o {output.vcf} -s {wildcards.sample} --summary-file {output.summary[0]} --stats-file {output.summary[1]} --blocks-file {output.summary[2]}
        '''
