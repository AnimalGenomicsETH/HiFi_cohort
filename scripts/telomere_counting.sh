for i in *.telo; do awk -v S=${i%.telo} '$1~/[0-9]/&&$2~/[0-9]/ {print S,"pubescent",$2}' ../logs/find_telomeres/sample-${i%.telo}*.err; done >> metadata.csv
for i in *.telo; do awk -v S=${i%.telo} '{print S,$4-$3,$2}' $i; done
