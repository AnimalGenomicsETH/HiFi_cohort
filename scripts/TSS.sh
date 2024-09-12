# from https://academic.oup.com/g3journal/article/13/8/jkad108/7175390#413067665

wget https://api.faang.org/files/trackhubs/BOVREG_CAGE_EUROFAANG/ARS-UCD1.2/tissues_TSS/TSS_testis_countAbove10_bb_minusY.Bigbed
bigBedToBed -tsv TSS_testis_countAbove10_bb_minusY.Bigbed /dev/stdout |\
cut -f -4 | sed 's/chr//g' | sed '1d' | sort -k1,1V -k2,2n > TSS_sites.bed

wget https://api.faang.org/files/trackhubs/BOVREG_CAGE_EUROFAANG/ARS-UCD1.2/tissues_TSS-Enhancers/BC_testis_countAbove10_bb_minusY.Bigbed
bigBedToBed -tsv BC_testis_countAbove10_bb_minusY.Bigbed /dev/stdout |\
cut -f -4 | sed 's/chr//g' | sed '1d' | sort -k1,1V -k2,2n > TSS_enhancer.bed 
