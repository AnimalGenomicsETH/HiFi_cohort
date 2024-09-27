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
        {{ echo -e "chrom\\tstart\\tend\\tcpg_label" ; bedtools makewindows -g {input.fai} -w 1000 | awk -v OFS='\\t' '{{print $0,$1"_"$2}}' ; }} > {output}
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
        mem_mb = lambda wildcards: 5000 if wildcards.mode == 'islands' else 85000,
        walltime = lambda wildcards: '30m' if wildcards.mode == 'islands' else '2h'
    shell:
        '''
        methbat build --input-collection <(echo -e "{params.collection}") --output-profile {output.profile}
        '''

rule methbat_compare:
    input:
        rules.methbat_build.output['profile']
    output:
        'methylation/cohort.{mode}.compare'
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '30m'
    shell:
        '''
        methbat compare --input-profile {input} --output-comparison {output} \
        --compare-category EPIDIDYMIS --baseline-category TESTIS \
        --min-zscore 3 --min-delta 0.2 --min-samples 10
        '''

rule compare_compares:
    input:
        expand(rules.methbat_compare.output,mode=('islands','TPM_testis.l0.01','TPM_testis.g5'))
    output:
        'methylation/comparisons.csv'
    localrule: True
    shell:
        '''
        echo "regions InsufficientData Uncategorized HyperASM HypoASM HyperMethylated HypoMethylated" > {output}
        for F in {input}
        do
          awk -v N=$(basename ${{F%.compare}}) 'NR>1 {{++L[$6]}} END {{print N,L["InsufficientData"],L["Uncategorized"],L["HyperASM"],L["HypoASM"],L["HyperMethylated"],L["HypoMethylated"]}}' $F
        done >> {output}
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

rule get_lowly_expressed_genes:
    input:
        TPM = lambda wildcards: '/cluster/work/pausch/HiFi_QTL/gene_expression/testis/UCD2.0_masked/117_samples/TPM_gene_testis_UNFILTERED.tsv' if wildcards.tissue == 'testis' else '/cluster/work/pausch/HiFi_QTL/gene_expression/testis/UCD2.0_masked/117_samples/TPM_gene_testis_UNFILTERED.tsv' #'/cluster/work/pausch/alex/Pop_HiFi/methylation/tissue_specific/Epidiymis_gene_TPM_unfiltered_UCD1.2.tsv'
    output:
        bed = 'methylation/tissue_specific/TPM_{tissue}.{logic}.list'
    params:
        logic = lambda wildcards: f'>{wildcards.logic[1:]}' if wildcards.logic[0] == 'g' else f'<{wildcards.logic[1:]}'
    localrule: True
    shell:
        '''
        awk -v OFS='\\t' 'NR>1 {{a=0; for (i=6;i<=NF;i++) {{a+=$i;N=(NF-6)}}; if(a/N {params.logic}) {{ print $4,$6,a/N,$1,$2,$3 }} }}' {input.TPM} > {output.bed}
        '''

rule get_gene_start_coordinates:
    input:
        gtf = 'GCF_002263795.3_ARS-UCD2.0_genomic.gtf.gz',
        genes = 'methylation/tissue_specific/{reason}.list',
        fai = config['reference'] + ".fai"
    output:
        bed = 'methylation/tissue_specific/{reason}.TSS.bed'
    localrule: True
    shell:
        '''
        zgrep -wf <(awk '{{print "\\""$1"\\""}}' {input.genes}) {input.gtf} | grep "start_codon" | cut -d' ' --output-delimiter=$'\\t' -f 1,4-5,7,10 | sort -u | awk '{{print $0"\\t{wildcards.reason}"}}' |\
        sort -k1,1V -k2,2n |\
        awk -v OFS='\\t' '$4=="+" {{ print $1,$2-1250,$3-750,$4;next }} {{print $1,$2+750,$3+1250,$4}}' > {output.bed} #adjust start coding to account for 500bp window around TSS, which is 1 Kb up/down of start_codon
        '''

# Data from https://academic.oup.com/g3journal/article/13/8/jkad108/7175390#413067665
rule prepare_TSS:
    output:
        multiext('methylation/TESTIS_TSS','.sites.bed','.enhancers.bed')
    localrule: True
    shell:
        '''
        TMPDIR=$(mktemp -d)
        wget -P $TMPDIR https://api.faang.org/files/trackhubs/BOVREG_CAGE_EUROFAANG/ARS-UCD1.2/tissues_TSS/TSS_testis_countAbove10_bb_minusY.Bigbed
        bigBedToBed -tsv $TMPDIR/TSS_testis_countAbove10_bb_minusY.Bigbed /dev/stdout |\
        cut -f -4 | sed 's/chr//g' | sed '1d' | sort -k1,1V -k2,2n > {output[0]}

        wget -P $TMPDIR https://api.faang.org/files/trackhubs/BOVREG_CAGE_EUROFAANG/ARS-UCD1.2/tissues_TSS-Enhancers/BC_testis_countAbove10_bb_minusY.Bigbed
        bigBedToBed -tsv $TMPDIR/BC_testis_countAbove10_bb_minusY.Bigbed /dev/stdout |\
        cut -f -4 | sed 's/chr//g' | sed '1d' | sort -k1,1V -k2,2n > {output[1]}
        '''

rule reformat_TSS:
    input:
        rules.get_gene_start_coordinates.output['bed']
    output:
        'methylation/CpG_{reason}.bed'
    localrule: True
    shell:
        '''
        {{ echo -e "chrom\\tstart\\tend\\tcpg_label" ; cat {input} ; }} > {output}
        '''

## CpG site analysis

#this will be a lot of data, so save encoding redundant information in the filename rather than per-line
rule convert_bed_to_csv:
    input:
        expand(rules.pb_CpG_tools.output['bed'],mapper='mm2',allow_missing=True)
    output:
        'methylation/sites/{sample}.{chromosome}.csv.gz'
    threads: 1
    resources:
        mem_mb = 1500,
        walltime = '30m'
    shell:
        '''
        {{ echo "position {wildcards.sample}" ; awk '$1=={wildcards.chromosome} {{ print $2,$4 }}' {input} ; }} | pigz -p 2 -c > {output}
        '''

rule find_all_reference_CpG_dinucleotides:
    input:
        reference = config['reference']
    output:
        'methylation/sites/reference.bed'
    localrule: True
    run:
        import regex
        import gzip
        CpG = regex.compile('CG')
        with open(output[0],'w') as fout, gzip.open(input.reference,'rt') as fin:
            offset = 0
            ends_in_C = False
            for line in (line.rstrip().upper() for line in fin):
                if line.startswith('>'):
                    chromosome = line[1:]
                    offset = 0
                    ends_in_C = False
                else:
                    if ends_in_C and line.startswith('G'):
                        fout.write(f'{chromosome}\t{offset-1}\t{offset+1}\n')
                    for hit in CpG.finditer(line):
                        fout.write(f'{chromosome}\t{offset+hit.start()}\t{offset+hit.end()}\n')
                    offset += len(line)
                    ends_in_C = line.endswith('C')

rule print_GT_matrix:
    input:
        vcf = multiext('mm2_DV/{chromosome}.Unrevised.vcf.gz','','.tbi'),
        CGs = rules.find_all_reference_CpG_dinucleotides.output
    output:
        'methylation/sites/{chromosome}.GTs.csv.gz'
    threads: 1
    resources:
        mem_mb = 2500,
        walltime = '1h'
    shell:
        '''
        {{ echo -n "position REF ALT " ; bcftools query -l {input.vcf[0]} | tr '\\n' ' ' | sed 's/.$/\\n/' ; bcftools query -f '%POS %REF %ALT[ %GT]' {input.vcf[0]} ; }} | pigz -p 2 -c > {output}
        '''
