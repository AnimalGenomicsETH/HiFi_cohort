awk '{print "SV",$2,$3}' /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt 
grep -Ff <(grep -Ff <(cut -f 1 /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt) final_set/eqtl_conditionals.txt | cut -d ' ' -f 8) /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt | awk '{print "eQTL",$2,$3}' 
grep -Ff <(grep -Ff <(cut -f 1 /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt) final_set/sqtl_conditionals.txt | cut -d ' ' -f 10) /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt | awk '{print "sQTL",$2,$3}'
grep -Ff <(bcftools +split-vep -f '%ID %IMPACT' refseq_vep.vcf.gz| grep "HIGH" | cut -d' ' -f 1) /cluster/work/pausch/HiFi_QTL/QTL/variants/117_samples/sv_length.txt | awk '{print "VEP",$2,$3}'
