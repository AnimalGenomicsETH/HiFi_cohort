rule fibertools_predict_m6a:
    input:
        'alignments/uBAM/{sample}/{cell}.5mC.bam'
        #raw PacBio ubam with kinetics
    output:
        'alignments/uBAM/{sample}/{cell}.m6a.bam' #new ubam with methylation but without kinetics
    threads: 16
    resources:
        mem_mb = 1500,
        walltime = '24h'
    shell:
        '''
        ft m6a -t {threads} {input} {output}
        '''

rule pb_CpG_tools:
    input:
        bam = rules.hiphase.output['bam']
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
        expand(rules.pb_CpG_tools.output['bed'],sample=samples,mapper='mm2')
    output:
        'methylation/merged.bed'
    shell:
        '''
        methbat profile --input-prefix {params.prefix} --output-region-profile {output}
        '''

rule methbat_build:
    input:
        collection = '',
        profiles = expand(rules.methbat_profile.output)
    output:
        ''
    shell:
        '''
        methbat build --input-collection {input.collection} --output-profile {output.profile}
        '''

rule methbat_compare:
    input:
        rules.methbat_build.output
    output:
        ''
    shell:
        '''
        methbat compare --input-profile {input} --output-comparison {output}
        '''
