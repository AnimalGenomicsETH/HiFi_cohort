#!/bin/bash

i=$1
bcftools query -f '%ID\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN' ${i}.hex.vcf.gz | awk '$4!="."&&$3!~/</ {if (length($2)>length($3)) {print ">"$1"_"$4"_"$5"\n"$2} else {print ">"$1"_"$4"_"$5"\n"$3}}' | seqtk seq -U - > ${i}.SVs.fa

RepeatMasker -pa 4 -e rmblast -lib /cluster/work/pausch/alex/REF_DATA/Repeat_libraries/BosTau9_repeat_library.fasta -q -no_is ${i}.SVs.fa

awk 'NR>2 {print $5,$7-$6,$10,$11}' ${i}.SVs.fa.out > ${i}.SVs.repeats

TRF ${i}.SVs.fa 2 5 7 80 10 50 2000 -h -ngs > ${i}.SVs.fa.trf


echo "SV length element class" > SVs.repeats
echo "SV period L score motif" > SVs.VNTRs
for i in {1..29}
do
  awk 'NR>3 {print $5,$7-$6,$10,$11}' ${i}.SVs.fa.out >> SVs.repeats
  awk '{if ($1~/@/) {S=substr($1,2); next} {print S,$3,$4,$8,$14}}' ${i}.SVs.fa.trf >> SVs.VNTRs
done
