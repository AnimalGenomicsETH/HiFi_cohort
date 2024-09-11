from pathlib import PurePath

rule pb_CpG_tools:
    input:
        bam = multiext('phasing/{sample}.mm2.phased.cram','','.cram')
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
        aligned_bam_to_cpg_scores --bam {input.bam[0]} --ref {config[reference]} --model {params.model} --output-prefix {params.prefix} --threads {threads}
        '''

rule prepare_CpGs:
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

rule methbat_profile:
    input:
        bed = expand(rules.pb_CpG_tools.output['bed'],mapper='mm2',allow_missing=True),
        CpGs = rules.prepare_CpGs.output
    output:
        profile = 'methylation/{sample}.BAT.profile',
        ASM = 'methylation/{sample}.BAT.ASM'
    params:
        prefix = lambda wildcards, input: PurePath(input['bed'][0]).with_suffix('').with_suffix('')
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '30m'
    shell:
        '''
        methbat profile --input-prefix {params.prefix} --input-regions {input.CpGs} --output-region-profile {output.profile} --output-asm-bed {output.ASM}
        '''

samples_M = [l.strip() for l in open('config/Braunvieh_testis_subset.txt')]

rule methbat_build:
    input:
        collection = expand(rules.methbat_profile.output['profile'],sample=samples_M)
    output:
        profile = 'methylation/cohort.BAT'
    params:
        collection = lambda wildcards, input: 'identifier\\tfilename\\tlabels\\n'+'\\n'.join(f'{S}\\t{P}\\tMALE' for S,P in zip(samples_M,input.collection))
    resources:
        mem_mb = 15000
    shell:
        '''
        methbat build --input-collection <(echo -e "{params.collection}") --output-profile {output.profile}
        '''

## re-profile using background as well as compare
rule methbat_compare:
    input:
        bed = expand(rules.pb_CpG_tools.output['bed'],mapper='mm2',allow_missing=True),
        profile = rules.methbat_profile.output.profile,
        background = rules.methbat_build.output
    output:
        profile = 'methylation/{sample}.BAT.cohort.profile',
        ASM = 'methylation/{sample}.BAT.cohort.ASM'
    params:
        prefix = lambda wildcards, input: PurePath(input['bed'][0]).with_suffix('').with_suffix('')
    resources:
        mem_mb = 5000,
        walltime = '30m'
    shell:
        '''
        methbat profile --input-prefix {params.prefix} --input-regions {input.background} --output-region-profile {output.profile} --output-asm-bed {output.ASM}
        #methbat compare --input-profile {input} --output-comparison {output}
        '''

rule methbat_gather:
    input:
        expand(rules.methbat_compare.output,sample=samples_M)
    output:
        'methylation/cohort.profiles.csv.gz'
    localrule: True
    shell:
        '''
        awk 'NR==1&&FNR==1 {{print "sample\\t"$0; next}} {{if (NR>1) {{split(FILENAME,a,"."); print a[1]"\\t"$0 }} }}' {input} | pigz -p 2 > {output}
        '''
