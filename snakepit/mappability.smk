from pathlib import PurePath

rule genmap_index:
    input:
        reference = config['reference']
    output:
        _dir = directory('mappability/index')
    threads: 1
    resources:
        mem_mb = 35000
    shell:
        '''
        genmap index -F {input.reference} -I {output._dir} -A divsufsort
        '''

#double check sanity that mappability is *improved* by longer k-mer
rule genmap_map:
    input:
        index = rules.genmap_index.output
    output:
        _dir = 'mappability/{kmer}_{error}.bedgraph'
    params:
        basename = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 1
    resources:
        mem_mb = 25000,
        walltime = '24h'
    shell:
        '''
        genmap map -K {wildcards.kmer} -E {wildcards.error} -I {input.index} -O {params.basename} -bg
        '''

rule bedtools_mappability:
    input:
        bed = expand(rules.genmap_map.output,error=(1,),kmer=(31,251)),
        fai = config['reference'] + ".fai"
    output:
        bed = expand('mappability/{kmer}_{error}.bed',error=(1,),kmer=(31,251)),
        regions = expand('mappability/{region}.bed',region=('easy','long','hard'))
    localrule: True
    shell:
        '''
        awk -v OFS='\\t' '$4>0.9 {{print $1,$2,$3}}' {input.bed[0]} | bedtools merge -i /dev/stdin > {output.bed[0]}
        awk -v OFS='\\t' '$4>0.9 {{print $1,$2,$3}}' {input.bed[1]} | bedtools merge -i /dev/stdin > {output.bed[1]}


        # get regions that were always okay
        bedtools intersect -a {output.bed[1]} -b {output.bed[0]} > {output.regions[0]}
        # get regions that are now more accessible
        bedtools subtract -a {output.bed[1]} -b {output.bed[0]} > {output.regions[1]}
        # get regions that are still hard
        bedtools complement -i {output.bed[1]} -g {input.fai} > {output.regions[2]}
        '''
