### Amova
this script performs a hierarchical analysis of variance as suggested during the review process.

original code from [https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html](https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html)

### Load packages
I installed them using conda
```r
library("vcfR")
library("adegenet")
library("poppr")
```

### Do the amova
```r
a<-read.vcfR("~/Desktop/bully_gbs/output_files/populations.snps.vcf")
data<-vcfR2genlight(a)
strata<-read.table("~/Desktop/bully_gbs/metadata/metadata_anova.txt",h=T,row.names="Sample")
strata(data)<-strata
amo<- poppr.amova(data,~Lake/depth)
```

```
$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                              Df     Sum Sq   Mean Sq
Between Lake                   1   	.590 8204.5902
Between depth Within Lake      8   4899.302  612.4128
Between samples Within depth  84  51393.022  611.8217
Within samples                94  57842.000  615.3404
Total                        187 122338.915  654.2188

$componentsofcovariance
                                                Sigma             %
Variations  Between Lake                  80.91369876  11.650190122
Variations  Between depth Within Lake      0.03209437   0.004621041
Variations  Between samples Within depth  -1.75936621  -0.253318673
Variations  Within samples               615.34042553  88.598507510
Total variations                         694.52685246 100.000000000

$statphi
                            Phi
Phi-samples-total  1.140149e-01
Phi-samples-depth -2.867374e-03
Phi-depth-Lake     5.230392e-05
Phi-Lake-total     1.165019e-01
```

My understanding is that there are more variance between lakes than between depths (almost none) but that there is considerable variance within samples that indicates the fact that the population structure is not dramatically strong.

let's look at significance:

```r
Amosignif   <- randtest(amo, nrepet = 1000)
```

```
                        Test          Obs    Std.Obs   Alter      Pvalue
1  Variations within samples 615.34042553 -4.0469194    less 0.000999001
2 Variations between samples  -1.75936621 -0.1262750 greater 0.544455544
3   Variations between depth   0.03209437  0.5114881 greater 0.306693307
4    Variations between Lake  80.91369876  5.4022616 greater 0.000999001
```

Only variation between lakes is significant.

