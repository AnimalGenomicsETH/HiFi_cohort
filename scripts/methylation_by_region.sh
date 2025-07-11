#!/bin/bash

#methbat signature --threads 6 --baseline-category OBV --compare-category BSW --input-collection cohort_information_with_breeds.TESTIS.files.csv --output-prefix BREED_SPECIFIC

FLANK=1000

echo "sample region coordinate methylation comparator"
while read -r C S E G
do
  region="$C:$S-$E"

  for I in $(awk 'NR>1 {print $1}' cohort_information_with_breeds.TESTIS.csv)
  do
    awk -v C=$C -v S=$S -v E=$E -v I=$I -v R=$region -v G=$G '$1==C&&$2>=S&&$3<E {print I,R,$2,$4,G}' $I.mm2.combined.bed
  done
done < <(awk 'NR>1&&($3-$2)>100' BREED_SPECIFIC.signature_regions.bed)
