## Comparison between bwa-mem2 and strobealign

### Accuracy

For 120 samples, we align HiFi with minimap2 and Illumina with bwa-mem2 and strobealign.
We call variants with DeepVariant using PACBIO or WGS models respectively, and then assess F1 with happy, taking mm2 as truth and bwa or strobe as query.

Overall, bwa is closer to mm2 than strobe, but the differences are fairly minor, although do escale to a 0.01 and 0.035 net difference in F1 for the X and Y respectively.
![image](https://github.com/AnimalGenomicsETH/HiFi_cohort/assets/29678761/88e59fc2-31ae-4c71-9afb-19804f9d94e8)

### Compute

There are some minor inconsistencies as some samples were run with slightly more cores to finish in time, and also measures the time/memory **including** the piping to `samtools sort` etc.
However, on average strobe finishes 3x faster and uses 40% less RAM.

![image](https://github.com/AnimalGenomicsETH/HiFi_cohort/assets/29678761/54c4886f-32c9-4d70-b15d-b7aeced07cb1)
