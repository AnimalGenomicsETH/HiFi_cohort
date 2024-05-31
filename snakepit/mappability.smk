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
        _dir = directory('mappability/{kmer}_{error}')
    threads: 1
    resources:
        mem_mb = 75000
    shell:
        '''
        genmap map -K {wildcards.kmer} -E {wildcards.error} -I {input.index} -O {output._dir} -bg
        '''

rule bedtools_mappability:
    input:
        bed = expand(rules.genmap_map.output,error=(1,),kmer=(31,251)),
        fai = config['reference'] + ".fai"
    output:
        'mappability/region_breakdown.bed'
    localrule: True
    shell:
        '''
        # get regions that were always okay
        bedtools intersect -A {input.bed[1]} -B {input.bed[0]} ...
        # get regions that are now more accessible
        bedtools subtract -A {input.bed[1]} -B {input.bed[0]} ...
        # get regions that are still hard
        bedtools complement -A {input.bed[1]} -g {input.fai} ...
        '''
