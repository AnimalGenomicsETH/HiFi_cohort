# HiFi sequencing of a Braunvieh population

Constructing a high quality and nearly homogeneous cohort of primarily Braunvieh cattle with DNA and RNA sequencing has proven to be a rich resource to explore trait associated phenotypes.

 - [Molecular quantitative trait loci in reproductive tissues impact male fertility in cattle](https://www.nature.com/articles/s41467-024-44935-7)
 - [Pangenome-genotyped structural variation improves molecular phenotype mapping in cattle](https://genome.cshlp.org/content/34/2/300.short)
 - [RNA-DNA differences in variant calls from cattle tissues result in erroneous eQTLs](https://link.springer.com/article/10.1186/s12864-024-10645-z)

Here, we extend the cohort with HiFi sequencing, to more thoroughly explore the role of structural variants on molecular QTL.

## Usage

Snakemake pipelines covering the major steps of HiFi alignment, small and structural variant calling, and association mapping are included here.
Some input (e.g., RNA alignments) are already assumed to be present.
The three major steps are detailed below

 - SV analysis
 - small variant comparison
 - association mapping

![rulegraph](rulegraph.svg]

n.b. These pipelines are designed to run on the Euler cluster, and may implicitly assume the majority of tools/environments are on path.

## Citation

> *Structural variants are enriched for molecular QTL in a cattle long read cohort*. Mapel, X.M., Leonard, A.S., Pausch, H. (**2025**).
