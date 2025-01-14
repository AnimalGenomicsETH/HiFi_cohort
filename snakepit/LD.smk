from pathlib import PurePath

rule separate_small_and_structural_variants:
    input:
        vcf = '/cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/filtered_merged/filtered.{chromosome}.vcf.gz'
    output:
        small = multiext('LD/small.{chromosome}.vcf.gz','','.csi'),
        structural = multiext('LD/structural.{chromosome}.vcf.gz','','.csi'),
        merged = multiext('LD/merged.{chromosome}.vcf.gz','','.csi')
    shell:
        '''
        bcftools view --threads {threads} -e 'ID~"Sniffles"' {input.vcf} |\
        bcftools annotate --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' --write-index -o {output.small[0]}
        bcftools view --threads {threads} -i 'ID~"Sniffles"' --write-index -o {output.structural[0]} {input.vcf}
        bcftools concat --threads {threads} {output.small[0]} {output.structural[0]} |\
        bcftools sort -T $TMPDIR --write-index -o {output.merged[0]} -
        '''

rule plink_tagging:
    input:
        vcf = rules.separate_small_and_structural_variants.output['merged'][0]
    output:
        'LD/SV.{chromosome}.r2_{r2}.tags.list'
    threads: 4
    resources:
        mem_mb = 5000,
        walltime = '1h'
    params:
        mem = lambda wildcards, threads, resources: threads*resources.mem_mb,
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    shell:
        '''
        plink --vcf {input.vcf} --show-tags <(bcftools query -i 'ID~"Sniffles"' -f '%ID' {input.vcf}) --tag-r2 0.{wildcards.r2} --tag-kb 1000 --threads {threads} --memory {params.mem} --chr-set 30 --vcf-half-call h --list-all --out {params.out}
        '''

rule summarise_tags:
    input:
        expand(rules.plink_tagging.output,chromosome=list(map(str,range(1,30)))+['X_het_missing','Y_HAP_het_missing','Y_PAR'],allow_missing=True)
    output:
        'LD/SV.r2_{r2}.tags.csv'
    localrule: True
    shell:
        '''
        echo "ID N_tags N_SV_tags Span_kb" > {output}
        awk '$1!="SNP" {{print $1,$4,gsub("Sniffles","",$8),$7}}' {input} >> {output}
        '''

