i=$1
bcftools query -f '%ID\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN' ${i}.hex.vcf.gz | awk '$4!="."&&$3!~/</ {if (length($2)>length($3)) {print ">"$1"_"$4"_"$5"\n"$2} else {print ">"$1"_"$4"_"$5"\n"$3}}' | seqtk seq -U - > ${i}.SVs.fa
RepeatMasker -pa 2 -e rmblast -lib /cluster/work/pausch/alex/REF_DATA/Repeat_libraries/BosTau9_repeat_library.fasta -s -no_is ${i}.SVs.fa
TRF ${i}.SVs.fa 2 5 7 80 10 50 2000 -h -ngs > ${i}.SVs.fa.trf
