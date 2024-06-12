
C=$1
S=$2
E=$3
W=$4

A=$((S-W))
R=$((E+W))

echo -e "sample\tbefore\tin\tafter"
for i in /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_SR/*cram
do
  echo -ne "$(basename ${i%.bwa.cram})\t"; samtools bedcov --reference /cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa  <(echo -e "$C\t$A\t$S\n$C\t$S\t$E\n$C\t$E\t$R") $i | awk '{print $4/($3-$2)}' |tr '\n' '\t' 
  echo
done
