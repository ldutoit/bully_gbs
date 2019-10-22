# SNPcalling.md
Data from bullies for Lake wakatipu or wanaka.

## 1. Understanding/Exploring the data

The data is single end two lanes. To understand the structure, the adaptor content and the barcodes I subset the two big *.gz* files to 250'000 reads.

```
# in source_files
# module load FastQC
zcat  SQ0353_CA62JANXX_s_1_fastq.txt.gz | head -n 1000000 > sample.fq
```

I now check the quality of the sequencing run using fastqc 

```bash
fastqc *fq
```

Output file is  [metadata/sample_fastqc.html](metadata/sample_fastqc.html) 


Single-end  100bp. Sequence quality is okay but 60% of sequences have adapters in them, it is a lot.

I trim them to 65bp. I need one set length for STACKS to work efficiently.

```bash
module load cutadapt
cutadapt --length  65  -a AGATCGGAAGAGC  -m 65  -o trimmed_samplt.fastq   sample.fq

fastqc trimmed_*
```


Output file is  [metadata/trimmed_samplt_fastqc.html](metadata/trimmed_samplt_fastqc.html) and it looks great!

I therefore trim the big dataset!
``` bash
cutadapt --length  65  -a AGATCGGAAGAGC  -m 65  -o lane1.fastq    SQ0353_CA62JANXX_s_1_fastq.txt.gz 

fastqc lane1*
```
Output file is [metadata/lane1_fastqc.html](metadata/lane1_fastqc.html) We retain **149'582'581** reads.

I then try to read demultiplexing on those after creating a borcode file.

```bash
##Barcode files from source_files_folder 
 cat ~/repos/scripts/bully_gbs/metadata/metadata_clean.txt  | awk -F " " '{print $8 "\t" $1}'| tail -n 95 > barcodes_lane1.txt

#from rootfolder

mkdir raw samples
cd raw
ln -s ../source_files/BullyGBS/lane1.fastq . 
cd ..


#module load Stacks
process_radtags  -p raw/ -o ./samples/ -b barcodes_lane1.txt -e pstI -r -c -q --inline-null
```
It seems to work ok:

```
149582581 total sequences
  4386035 barcode not found drops (2.9%)
   135493 low quality read drops (0.1%)
   411357 RAD cutsite not found drops (0.3%)
144649696 retained reads (96.7%)
```


### Run denovo_map.pl full

This wrapper does all the matching between samples. We'll run it with the default parameters. Both lanes had the same individuals so we can pick any barcode files. We just make sure we exclude the negative control.

```bash
#from root folder create popmap
cut -f 2 barcodes_lane1.txt | awk '{print $0,"\tsinglepop"}'  - |  grep -v NEG > popmap_allNONEG.txt
mkdir output_M2
echo '#!/bin/sh' > cleanM2.sh
echo "denovo_map.pl --samples samples/ --popmap popmap_allNONEG.txt  -o output_M2  -M 2 -n 2 -m 3 -T 16"  >> cleanM2.sh
sbatch -A uoo00116 -t 1-00:00:00 --partition=long -J cleanM2 -c 16 --mem=64G cleanM2.sh # specific to the mahuika cluster and ludovic.dutoit

```

We then run populations excluding no samples and only the SNPs with more than 65% heterozygosity as likely collapsed paralogs and remove any SNP not covered in at least 80% of individuals.

```bash
 populations -P output_M2/ -M popmap_allNONEG.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -r 0.8  --write-single-snp # then filter it without
 ```


*That outputted 17314 SNPs* It is quite good, it means we have a lot of well-covered SNPs and therefore a lot of space for filtering.

**Exploration**

I then use a little bit of R code to investigate how much SNPs we have across how many individuals? We

```r
library("VariantAnnotation")
data<-readVcf("populations.snps.vcf")


numbermissing<-function(x){
	return(length(grep("\\./\\.",x)))
}

countsofmissing<-apply(geno(data)$GT,2,numbermissing)
nonmissing<-dim( geno(data)$GT)[1]-countsofmissing
### get list of removals ad then re run as white list with  proportion of missing look
quantile(nonmissing,seq(0,1,0.02))
      0%       2%       4%       6%       8%      10%      12%      14%
15539.00 16050.28 16124.44 16288.16 16346.76 16396.20 16441.32 16541.56
     16%      18%      20%      22%      24%      26%      28%      30%
16579.52 16605.80 16662.40 16670.68 16714.00 16718.44 16742.28 16747.40
     32%      34%      36%      38%      40%      42%      44%      46%
16758.00 16759.92 16768.20 16787.28 16797.60 16807.88 16821.88 16828.00
     48%      50%      52%      54%      56%      58%      60%      62%
16829.60 16852.00 16860.76 16887.56 16896.28 16899.60 16905.00 16912.40
     64%      66%      68%      70%      72%      74%      76%      78%
16916.00 16923.24 16930.84 16934.00 16947.40 16957.24 16961.44 16970.28
     80%      82%      84%      86%      88%      90%      92%      94%
16974.80 16979.16 16982.92 16987.68 16992.16 16996.00 17002.84 17009.36
     96%      98%     100%
17012.48 17017.12 17050.00
```

That actually looks great, with all individuals having at leas 15539 SNPs!

I am therefore happy with this SNP filtering that outputted 17314 SNPs in total.


## No missing data

to try to maximise PCA and faststructure resolution, I create a version without missing data.

```
vcftools --vcf populations.snps.vcf --max-missing 1.0 --recode 
mv out.recode.vcf populations.snps.NOMISSING.vcf
```

That is the final version of the SNP calling and the vcf as well as some other format go into [output_files/](output_files)