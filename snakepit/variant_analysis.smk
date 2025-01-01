rule plink2_LD:
    input:
        vcf = 'filtered.{chromosome}.vcf.gz'
    output:
        vcf = '{chromosome}.hex.vcf.gz'
    params:
        prefix = '{chromosome}_SVs'
    shell:
        '''
        bcftools annotate --set-id "%VKX" -o {output.vcf} -W=tbi {input.vcf}
        zgrep "SVTYPE=" {output.vcf} | cut -f 3 > $TMPDIR/SV_list
        plink2 --vcf {output.vcf} --out {params.prefix} --threads 2  --vcf-half-call missing --ld-snp-list $TMPDIR/SV_list --r2-unphased zs --ld-window-r2 0 --ld-window 999999999 --cow
        '''

rule estimate_DUP_genotypes:
    input:
        fai = config['reference']+'.fai',
        vcf = '/cluster/work/pausch/HiFi_QTL/final_set/final_filtered.all.vcf.gz'
    output:
        csv = 'duplications.csv'
    params:
        cram_path = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_HiFi/{}.mm2.cram'
    shell:
        '''
        awk -v OFS='\\t' 'NR<32 {{print $1,1,$2}}' {input.fai} > $TMPDIR/regions.bed
        bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE' {input.vcf} |\
        grep "DUP" |\
        cut -f -3 >> $TMPDIR/regions.bed

        echo "sample read chromosome start end coverage count" | sed 's/ /\\t/g' > {output.csv}

        parallel --line-buffer -j $1 "samtools bedcov --reference {config[reference]} -c $TMPDIR/regions.bed {params.cram_path} |\
        awk -v OFS='\\t' -v S={} '{print S,\\"HiFi\\",\$0}'; samtools bedcov --reference {config[decompress_reference]} -c $TMPDIR/regions.bed {params.cram_path} |\
        awk -v OFS='\t' -v S={} '{print S,\"SR\",\$0}'" ::: $(bcftools query -l {input.vcf}) >> {output.csv}

        echo -e "chrom\\tpos\\t$(bcftools query -l {input.vcf} | tr '\\n' '\\t')" | sed 's/\\t$//' > {output.GTs}
        
        bcftools query -i 'INFO/SVTYPE=="DUP"' -f '%CHROM\\t%POS[\\t%GT]' {input.vcf} >> {output.GTs}
        '''

#rule SV_vep?
#awk '{print "SV",$2,$3}' /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt 
#grep -Ff <(grep -Ff <(cut -f 1 /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt) final_set/eqtl_conditionals.txt | cut -d ' ' -f 8) /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt | awk '{print "eQTL",$2,$3}' 
#grep -Ff <(grep -Ff <(cut -f 1 /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt) final_set/sqtl_conditionals.txt | cut -d ' ' -f 10) /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt | awk '{print "sQTL",$2,$3}'
#grep -Ff <(bcftools +split-vep -f '%ID %IMPACT' refseq_vep.vcf.gz| grep "HIGH" | cut -d' ' -f 1) /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt | awk '{print "VEP",$2,$3}'
#bcftools +split-vep -f '%ID %SYMBOL|%IMPACT' refseq_vep.vcf.gz








