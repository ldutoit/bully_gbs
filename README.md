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

[structureexploration.md](structureexploration.md) makes a principal component analysis of the data, visualised by lake, exact depth or deep vs shallow. The data clearly cluster by lake, not so much by depth. I went on to do a faststructure analysis varying K between 2 to 5 with within lakes dataset or from K=2 to K=10 when grouping lakers together.  There is no structure whatsoever beyond lake.

NOTE: two individuals where clustering on their own at the first PCA andare actually the same individuals, they have been treated as such bringing the number of samples to 94

## Output files

All output_files are in [output_files](output_files). Note that samples metadata and barcodes information is in [metadata](metadata). Most of the figures are visualised directly from within [structurexploration.md](structurexploration.md), the pca score as stored in [output_files/scores_pca.txt] and the faststructureoutput are in [faststructure/](faststructure)

