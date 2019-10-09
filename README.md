# bully_gbs

## Description
SNP calling and structure + PCA analysis for bullies in two different lakes at different depths.

## Key Players

Travis Ingram

## Physical location of the data

contact Travis ngram for access.

## SNP Calling
The SNP calling is done SNPcalling.md with Stacks/2.41. 

It will be done on some variation of the [github.com/ldutoit/stoneflies_8pops](github.com/ldutoit/stoneflies_8pops) model (private, contact ldutoit for access)


*In construction*

Briefly, we started with x samples and a total of x paired/single end reads on x lanes of sequencing reads  not sure as of now on which instrument as of now). After checking the quality using FastQC v0.11.7, we removed adapter contamination and shortened reads to a common length of xbp using cutadaptv2.3. X reads remained. We demultiplexed those using the process_radtags of Stacks/v2.41. We then extracted SNPs de-novo using the denovo_map.pl wrapper of Stacks. We removed loci genotyped for less than 70% of individuals across all populations and kept one SNP per locus for a total of 2974 SNPs.

NOTE: X individuals (i.e. ... ) were removed for very poor sequencing (<200SNPs), and two extra because they were clear outliers in the tree ( potentially a different species leaving us with X individuals)

## Population structure analysis

*In construction*

## Output files

*In construction*
