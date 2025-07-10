
echo -e "mapper\tsample\ttime\tmemory"
for M in bwa strobe
do
  for S in $(cat sample_IDs.txt)
  do
    job_log=$(ls -Art logs/short_read_align/sample-${S}.mapper-${M}-202310*.err | tail -n 1)
    job_ID=$(grep -oP "scratch/tmp.\K\d+" ${job_log})
    myjobs -j ${job_ID} | awk -v M=$M -v S=$S 'BEGIN {printf "%s\t%s\t",M,S} /Total (CPU time|resident memory)/ {printf "%s\t",$5} END {print}' 
  done
done
