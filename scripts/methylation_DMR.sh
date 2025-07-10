#!/bin/bash

CHR=$1
START=$2
END=$3

echo "sample coordinate methylation"

for S in $(awk 'NR>1 {print $1}' cohort_information_with_breeds.TESTIS.csv)
do
  awk -v C=$CHR -v S=$START -v E=$END -v s=$S '$1==C&&$2>=S&&$3<E {print s,$2,$4}' $S.mm2.combined.bed
done
