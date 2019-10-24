# bully_gbs

## Description
SNP calling and structure + PCA analysis for bullies in two different lakes at different depths.

## Key Players

Travis Ingram

## Physical location of the data

contact Travis Ingram for access.

## SNP Calling
The SNP calling is done [SNPcalling.md](SNPcalling.md) with Stacks/2.41. 

 Done on some variation of the [github.com/ldutoit/stoneflies_8pops](github.com/ldutoit/stoneflies_8pops) model (private, contact ldutoit for access)


## Population structure analysis

[structurexploration.md](structurexploration.md) makes a principal component analysis of the data, visualised by lake, exact depth or deep vs shallow. The data clearly cluster by lake, not so much by depth. I went on to do a faststructure analysis varying K between 2 to 5 with within lakes dataset or from K=2 to K=10 when grouping lakers together.  There is no structure whatsoever beyond lake.

NOTE: two individuals where clustering on their own at the first PCA and have been excluded from further PCA and Structure analyses.

## Output files

All output_files are in [output_files](output_files). Note that samples metadata and barcodes information is in [metadata](metadata). Most of the figures are visualised directly from within [structurexploration.md](structurexploration.md)

## Methods in plain english


95 samples were sequenced for single-end 100bp reads on one lane of XXXXX. After checking the reads quality using FastQC v0.11.7, we removed adapter contamination and shortened reads to a common length of 65 bp using cutadaptv2.3. The remaining 149'582'581  reads were demultiplexed  using the process_radtags function  of Stacks/v2.41 (CIT). We then extracted SNPs de-novo using the denovo_map.pl wrapper of Stacks. Finally, we removed loci genotyped for less than 100% of individuals and kept one variant per locus for a total of 9138 SNPs. 
To investigate the relationship between lake depth and genetic structure,  Principal component analyses were performed using the pcadapt (CIT) package in R (CIT) and vcftools (CIT). Population structure between and within lakes was further investigates using fastStructure (CIT) and PGDSpider v.XXX for file conversion. K was varied between 2 to 5 representing the different depth within lake and k=2 to K=10 representing each possible combination of depth and lake for the analysis incorporating both lake.

Note: Add citations and software versions. Make it clear in the results that 2 inds (WK16-48 WK16-47) were removed and decide at which stage to do it.

