rule get_bad_sites:
    input:
        expand(rules.happy.output['others'][0],chromosome=main_regions,filtering=('filtered',),read1=('mm2',),read2=('bwa',),allow_missing=True)
    output:
        bed = 'eVariant_overlap/{sample}.bed'
    shell:
        '''
        bcftools query -e 'FORMAT/BD=="TP"' -f '%CHROM\\t%POS\\t%INFO/Regions\\n' -v <(ls {input}) | sort -u | awk -v OFS='\\t' '{{print $1,$2,$2+1,$3}}' | sort -k1,1n -k2,2n > {output.bed}
        '''

rule gather_bad_sites:
    input:
        expand(rules.get_bad_sites.output,sample=samples)
    output:
        'eVariant_overlap/bad_sites.bed'
    shell:
        '''
        cat {input} | sort -k1,1 -k2,2n | bedtools merge -i /dev/stdin -d -1 -c 4,4 -o first,count > {output}
        '''

rule overlap_bad_sites:
    input:
        QTL = expand('/cluster/work/pausch/HiFi_QTL/QTL/eQTL/117_samples/testing_HiFi_illumina/Testis_{read}_filtered/conditionals.{chromosome}.01.05.txt.gz',chromosome=range(1,30),allow_missing=True),
        bed = rules.gather_bad_sites.output
    output:
        'eVariant_overlap/bad_sites.{read}.csv'
    shell:
        '''
        zcat {input.QTL} |\
        awk -v OFS='\\t' '{{print $9,$10,$10,$19}}' |\
        bedtools intersect -wo -a /dev/stdin -b {input.bed} > {output}
        '''

