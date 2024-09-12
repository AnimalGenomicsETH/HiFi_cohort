from pathlib import PurePath

rule pb_CpG_tools:
    input:
        bam = multiext('/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_HiFi_phased/{sample}.mm2.phased.cram','','.crai'),
        reference = config['reference']
    output:
        bed = "methylation/{sample}.{mapper}.combined.bed",
        BigWig = "methylation/{sample}.{mapper}.combined.bw"
    threads: 4
    resources:
        mem_mb = 5000,
        walltime = '4h'
    params:
        prefix = lambda wildcards, output: PurePath(output['bed']).with_suffix('').with_suffix(''),
        model = '/cluster/work/pausch/alex/software/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite'
    shell:
        '''
        aligned_bam_to_cpg_scores --bam {input.bam[0]} --ref {input.reference} --model {params.model} --output-prefix {params.prefix} --threads {threads}
        '''

rule prepare_CpG_islands:
    output:
        'methylation/CpG_islands.bed'
    localrule: True
    shell:
        '''
        TMPDIR=$(mktemp -d)
        wget -P $TMPDIR https://hgdownload.soe.ucsc.edu/hubs/GCF/002/263/795/GCF_002263795.3/bbi/GCF_002263795.3_ARS-UCD2.0.cpgIslandExtUnmasked.bb
        wget -P $TMPDIR https://hgdownload.soe.ucsc.edu/hubs/GCF/002/263/795/GCF_002263795.3/GCF_002263795.3.chromAlias.txt

        bigBedToBed -tsv $TMPDIR/GCF_002263795.3_ARS-UCD2.0.cpgIslandExtUnmasked.bb /dev/stdout |\
        cut -f -4 | sed 's/ /_/g' | sed '1cchrom\tstart\tend\tcpg_label' |\
        awk -v OFS='\\t' 'NR==FNR{{a[$1]=$2;b[$1]=$3;next}}$1 in a{{$1=($1~/NC/?a[$1]:b[$1])}}1' $TMPDIR/GCF_002263795.3.chromAlias.txt - > {output}
        '''

rule prepare_CpG_windows:
    input:
        fai = config['reference'] + ".fai"
    output:
        'methylation/CpG_windows.bed',
    localrule: True
    shell:
        '''
        bedtools makewindows -g <({input.fai}) -w 1000 > {output}
        '''

rule methbat_profile:
    input:
        bed = expand(rules.pb_CpG_tools.output['bed'],mapper='mm2',allow_missing=True),
        regions = lambda wildcards:  'methylation/CpG_{mode}.bed' if wildcards._group == 'individual' else 'methylation/cohort.{mode}.BAT'
    output:
        BAT = multiext('methylation/{sample}.{mode}.{_group}.BAT','.profile','.ASM')
    params:
        prefix = lambda wildcards, input: PurePath(input['bed'][0]).with_suffix('').with_suffix('')
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '1h'
    shell:
        '''
        methbat profile --input-prefix {params.prefix} --input-regions {input.regions} --output-region-profile {output.BAT[0]} --output-asm-bed {output.BAT[1]}
        '''

def build_information_collection(samples,profiles,metadata):
    labels = {line.split()[0]:line.split()[2] for n,line in enumerate(open(metadata)) if n>0}
    return "identifier\\tfilename\\tlabels\\n" + '\\n'.join(f'{S}\\t{P}\\t{labels[S]}' for S,P in zip(samples,profiles))

rule methbat_build:
    input:
        collection = expand(rules.methbat_profile.output['BAT'][0],_group='individual',sample=samples,allow_missing=True),
        metadata = 'methylation/cohort_information.tsv'
    output:
        profile = 'methylation/cohort.{mode}.BAT'
    params:
        collection = lambda wildcards, input: build_information_collection(samples,input.collection,input.metadata)
    resources:
        mem_mb = 15000
    shell:
        '''
        methbat build --input-collection <(echo -e "{params.collection}") --output-profile {output.profile}
        '''

rule methbat_gather:
    input:
        expand(rules.methbat_profile.output['BAT'][0],_group='cohort',sample=samples,allow_missing=True)
    output:
        'methylation/{mode}.cohort.profiles.csv.gz'
    localrule: True
    shell:
        '''
        awk 'NR==1&&FNR==1 {{print "sample\\t"$0; next}} {{if (NR>1) {{split(FILENAME,a,"."); print a[1]"\\t"$0 }} }}' {input} | pigz -p 2 > {output}
        '''

rule tissue_specific_regions:
    input:
        gtf = 'GCF_002263795.3_ARS-UCD2.0_genomic.gtf',
        testis_specific_genes = 'methylation/tissue_specific/testis.list',
        epididymis_specific_genes = 'methylation/tissue_specific/epididymis.list'
    output:
        'methylation/tissue_specific/testis.bed'
    shell:
        '''
        grep -wf <(awk '{{print "\\""$1"\\""}}' {input.testis_specific_genes}) {input.gtf} | grep "start_codon" | cut -d' ' --output-delimiter=$'\\t' -f 1,4-5,7,10 | sort -u
        grep -wf <(awk '{{print "\\""$1"\\""}}' {input.epididymis_specific_genes}) {input.gtf} | grep "start_codon" | cut -d' ' --output-delimiter=$'\\t' -f 1,4-5,7,10 | sort -u
        '''
