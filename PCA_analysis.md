Principle component analysis (PCA) and Differential gene expression
-------------------------------------------------------------------

DEseq analysis based on the negative binomial distribution. Estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.

In this analysis, we apply Deseq2 into datasets collected from 58 RNA-Sequencing libraries from our lab and public data from NCBI to study co-expression network in light signaling in Arabidopsis First of all, we use PCA analysis for dimentional reduction analysis to analyze the quality of 58 RNA-seq samples

Please cite: Pham, V.N., Xu X, Huq E. Molecular bases for the constitutive photomorphogenic phenotypes in Arabidopsis. Development 2018:dev.169870 doi: 10.1242/dev.169870.

``` r
library(tidyverse)
theme_set(theme_bw(base_size=12)) # set default ggplot2 theme
```

First, to read our light condition data file "treatment" contained all the information about sample, condition, dataset and light conditions

``` r
samples <- read.csv(file.path("treatment.csv"), header=TRUE)
head(samples)
```

    ##       sample  condition dataset whole.seedlings Only.cotyledon
    ## 1 SRR1292205 Red 5 days       1               1              0
    ## 2 SRR1292206 Red 5 days       1               1              0
    ## 3 SRR1292207 Red 5 days       1               1              0
    ## 4 SRR1307153     FR 3h        2               1              0
    ## 5 SRR1421920 Blue 4days       3               1              0
    ## 6 SRR1421921 Blue 4days       3               1              0
    ##   Only.hypocotyl stage..days. Red.5.days Far.red.3h Blue.4days Low.Blue.1h
    ## 1              0            5          1          0          0           0
    ## 2              0            5          1          0          0           0
    ## 3              0            5          1          0          0           0
    ## 4              0            4          0          1          0           0
    ## 5              0            4          0          0          1           0
    ## 6              0            4          0          0          1           0
    ##   White.Light.1h Low.Blue.6h White.Light.6h Low.Blue.24h White.Light.24h
    ## 1              0           0              0            0               0
    ## 2              0           0              0            0               0
    ## 3              0           0              0            0               0
    ## 4              0           0              0            0               0
    ## 5              0           0              0            0               0
    ## 6              0           0              0            0               0
    ##   White.Light.18.days High.Light.20s High.Light.60s Low.Blue.16h
    ## 1                   0              0              0            0
    ## 2                   0              0              0            0
    ## 3                   0              0              0            0
    ## 4                   0              0              0            0
    ## 5                   0              0              0            0
    ## 6                   0              0              0            0
    ##   White.Light.16h Cotyledon.Dark Hypocotyl.Dark Cotyledon.Light.1h
    ## 1               0              0              0                  0
    ## 2               0              0              0                  0
    ## 3               0              0              0                  0
    ## 4               0              0              0                  0
    ## 5               0              0              0                  0
    ## 6               0              0              0                  0
    ##   Hypocotyl.Light.1h Cotyledon.Light.6h Hypocotyl.Light.6h Dark.4.days
    ## 1                  0                  0                  0           0
    ## 2                  0                  0                  0           0
    ## 3                  0                  0                  0           0
    ## 4                  0                  0                  0           0
    ## 5                  0                  0                  0           0
    ## 6                  0                  0                  0           0
    ##   Red.3h Red.2.days White.Light.30.days Red.1h
    ## 1      0          0                   0      0
    ## 2      0          0                   0      0
    ## 3      0          0                   0      0
    ## 4      0          0                   0      0
    ## 5      0          0                   0      0
    ## 6      0          0                   0      0

Import transcript abundance estimates from Kallisto transctip abundance files

``` r
library("tximport")
library("readr")
library("tximportData")
rownames(samples) <- samples$sample
files <- file.path("kallisto", samples$sample, "abundance.h5")
names(files) <- samples$sample
names(files)
```

    ##  [1] "SRR1292205"   "SRR1292206"   "SRR1292207"   "SRR1307153"  
    ##  [5] "SRR1421920"   "SRR1421921"   "SRR1523303"   "SRR1523304"  
    ##  [9] "SRR1523305"   "SRR1523306"   "SRR1523311"   "SRR1523312"  
    ## [13] "SRR1523313"   "SRR1523314"   "SRR1523319"   "SRR1523320"  
    ## [17] "SRR1523321"   "SRR1523322"   "SRR1560610"   "SRR1560611"  
    ## [21] "SRR1560612"   "SRR1560613"   "SRR1560614"   "SRR1560615"  
    ## [25] "SRR1560616"   "SRR1560617"   "SRR1560618"   "SRR3046818"  
    ## [29] "SRR3046819"   "SRR3046820"   "SRR3046821"   "SRR3290779"  
    ## [33] "SRR3290780"   "SRR3290781"   "SRR3290782"   "SRR3290784"  
    ## [37] "SRR3290785"   "SRR3290786"   "SRR3290788"   "SRR3290791"  
    ## [41] "SRR3290792"   "SRR3290794"   "SRR3290795"   "SRR3881040_1"
    ## [45] "SRR3881040_2" "SRR3881040_3" "SRR4046133"   "SRR4046134"  
    ## [49] "SRR4046135"   "SRR5020734"   "SRR5020735"   "SRR5020736"  
    ## [53] "SRR6215007"   "SRR6215008"   "SRR6215009"   "SRR7521802"  
    ## [57] "SRR7521803"   "SRR7521804"

``` r
conditions <-  c(rep("Red 5 days",3), rep("FR 3h ",1), rep("Blue 4days",2), rep("Low Blue 1h",2), rep("White Light 1h", 2), 
                 rep("Low Blue 6h",2), rep("White Light 6h", 2), rep("Low Blue 24h", 2), rep("White Light 24h", 2), rep("White Light 18 days",3), 
                 rep("High Light 20s",3), rep("High Light 60s", 3), rep("Low Blue 16h", 2), rep("White Light 16h", 2), rep("Cotyledon Dark", 3), 
                 rep("Hypocotyl Dark", 2), rep("Cotyledon Light 1h", 2), rep("Hypocotyl Light 1h", 2), rep("Cotyledon Light 6h", 1), rep("Hypocotyl Light 6h", 2),
                 rep("Dark 4 days", 3), rep("Red 3h,", 3), rep("Red 2 days", 3), rep("White Light 30 days", 3), rep("Red 1h",3))
                 
sampleTable <- data.frame( condition = conditions)
sampleTable
```

    ##              condition
    ## 1           Red 5 days
    ## 2           Red 5 days
    ## 3           Red 5 days
    ## 4               FR 3h 
    ## 5           Blue 4days
    ## 6           Blue 4days
    ## 7          Low Blue 1h
    ## 8          Low Blue 1h
    ## 9       White Light 1h
    ## 10      White Light 1h
    ## 11         Low Blue 6h
    ## 12         Low Blue 6h
    ## 13      White Light 6h
    ## 14      White Light 6h
    ## 15        Low Blue 24h
    ## 16        Low Blue 24h
    ## 17     White Light 24h
    ## 18     White Light 24h
    ## 19 White Light 18 days
    ## 20 White Light 18 days
    ## 21 White Light 18 days
    ## 22      High Light 20s
    ## 23      High Light 20s
    ## 24      High Light 20s
    ## 25      High Light 60s
    ## 26      High Light 60s
    ## 27      High Light 60s
    ## 28        Low Blue 16h
    ## 29        Low Blue 16h
    ## 30     White Light 16h
    ## 31     White Light 16h
    ## 32      Cotyledon Dark
    ## 33      Cotyledon Dark
    ## 34      Cotyledon Dark
    ## 35      Hypocotyl Dark
    ## 36      Hypocotyl Dark
    ## 37  Cotyledon Light 1h
    ## 38  Cotyledon Light 1h
    ## 39  Hypocotyl Light 1h
    ## 40  Hypocotyl Light 1h
    ## 41  Cotyledon Light 6h
    ## 42  Hypocotyl Light 6h
    ## 43  Hypocotyl Light 6h
    ## 44         Dark 4 days
    ## 45         Dark 4 days
    ## 46         Dark 4 days
    ## 47             Red 3h,
    ## 48             Red 3h,
    ## 49             Red 3h,
    ## 50          Red 2 days
    ## 51          Red 2 days
    ## 52          Red 2 days
    ## 53 White Light 30 days
    ## 54 White Light 30 days
    ## 55 White Light 30 days
    ## 56              Red 1h
    ## 57              Red 1h
    ## 58              Red 1h

Next we read in a table that links transcripts to genes for this dataset.

``` r
tx2gene <- read.csv("transcript_gene_mapping_full.csv")
```

Import the necessary quantification data for DESeq2 using the tximport function.

``` r
txi <- tximport(files, type="kallisto", tx2gene=tx2gene)
```

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

Then we can construct a DESeqDataSet from the txi object and sample information in samples.

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colMeans, colnames, colSums, dirname, do.call, duplicated,
    ##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
    ##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
    ##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
    ##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## Loading required package: BiocParallel

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply

``` r
deseq <- DESeqDataSetFromTximport(txi,
                                  colData = samples,
                                  design = ~ condition)
```

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not an warning or error]

    ## using counts and average transcript lengths from tximport

The standard differential expression analysis steps are wrapped into a single function, DESeq.

``` r
deseq
```

    ## class: DESeqDataSet 
    ## dim: 27655 58 
    ## metadata(1): version
    ## assays(2): counts avgTxLength
    ## rownames(27655): AT1G01010 AT1G01020 ... ATMG01400 ATMG01410
    ## rowData names(0):
    ## colnames(58): SRR1292205 SRR1292206 ... SRR7521803 SRR7521804
    ## colData names(32): sample condition ... White.Light.30.days Red.1h

``` r
keep <- rowSums(counts(deseq)) >= 10
deseq <- deseq[keep,]
deseq <- DESeq(deseq)
```

    ## estimating size factors

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not an warning or error]

    ## using 'avgTxLength' from assays(dds), correcting for library size

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not an warning or error]

    ## final dispersion estimates

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not an warning or error]

    ## fitting model and testing

    ## 1 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest

Variance stabilizing transformation

``` r
vsdB <- varianceStabilizingTransformation(deseq)
```

Data quality assessment by sample clustering and visualization
--------------------------------------------------------------

### Heatmap of the count matrix

To explore a count matrix, it is often instructive to look at it as a heatmap. Below we show how to produce such a heatmap for variant transfromation data from deseq

``` r
library("pheatmap")
select <- order(rowMeans(counts(deseq,normalized=TRUE)),
                decreasing=TRUE)[1:200]


df <- as.data.frame(colData(deseq)[,c("condition", "dataset")])

pheatmap(assay(vsdB)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

<img src="PCA_analysis_files/figure-markdown_github/unnamed-chunk-10-1.png" style="display: block; margin: auto;" /> Heatmap of the sample-to-sample distances.

Another use of the transformed data is sample clustering. Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.

``` r
library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(assay(vsdB)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsdB$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

<img src="PCA_analysis_files/figure-markdown_github/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

PCA analysis
------------

Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

``` r
pca<-  plotPCA(vsdB,intgroup = c ("condition"), returnData = TRUE) 
```

``` r
percentVar <- round(100 * attr(pca, "percentVar"))
```

``` r
cbPalette <- c("red", "blue", "violetred", "black", "pink", "purple", "cyan", "green", "yellow", "orange","magenta", "chocolate", "tan", 
               "turquoise", "midnightblue", "greenyellow", "darkcyan", "coral4", "lightseagreen", "maroon", "khaki", "salmon", "tomato1", "blue4", "deeppink4")
```

``` r
p <-  ggplot(pca, aes(PC1, PC2, color= condition, label = condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  scale_colour_manual(values=cbPalette)

p
```

<img src="PCA_analysis_files/figure-markdown_github/unnamed-chunk-15-1.png" style="display: block; margin: auto;" /> Instead of using condition to separate the samples, we wanted to see if there is any outliers in datasets or whether there is any bacth effect, so we group vcdB into group by dataset and plot the PCA.

``` r
pca<-  plotPCA(vsdB,intgroup = c ("dataset"), returnData = TRUE) 

cbPalette <- c("red", "blue", "violetred", "black", "pink", "purple", "cyan", "green", "yellow", "orange")
# PCA plot using ggplot  to customize PCA plot 

p <-  ggplot(pca, aes(PC1, PC2, color= factor(dataset), label = factor(dataset))) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  scale_colour_manual(values=cbPalette)

p
```

<img src="PCA_analysis_files/figure-markdown_github/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

Another way to perform PCA that we use prcomp function

``` r
pca <- prcomp(t(assay(vsdB)))
```

``` r
dataGG = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], 
                    PC3 = pca$x[,3], PC4 = pca$x[,4], 
                    sampleNO = vsdB$sample,
                    condition = vsdB$condition,
                    dataset = vsdB$dataset)
```

``` r
library(ggrepel)
```

    ## Warning: package 'ggrepel' was built under R version 3.5.2

``` r
percentVar <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
p <-  ggplot(dataGG, aes(PC1, PC2, color= factor(dataset), label = condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
p <- p + geom_text_repel(box.padding = 0.1) + scale_radius(range = c(3,6))
p
```

<img src="PCA_analysis_files/figure-markdown_github/unnamed-chunk-19-1.png" style="display: block; margin: auto;" /> The samples from dataset 6 (hypocotyl Dark, Light and cotyledon Dark, Light) spread out over PC1, but does cluster tightly along PC2

Conesa et al 2016 (<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8>)

Reproducibility among technical replicates should be generally high (Spearman R2 &gt; 0.9) \[1\], but no clear standard exists for biological replicates, as this depends on the heterogeneity of the experimental system."

``` r
p2 <-  ggplot(dataGG, aes(PC2, PC3, color= factor(dataset), label = condition)) +
  geom_point(size=4) +
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) + 
  coord_fixed()
p2 <- p2 + geom_text_repel(box.padding = 0.1) + scale_radius(range = c(3,6))
p2
```

<img src="PCA_analysis_files/figure-markdown_github/unnamed-chunk-20-1.png" style="display: block; margin: auto;" /> So we can see samples are separated by dataset better on PC2 and PC3. So it might be

``` r
percentVar <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
barplot(percentVar, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
```

<img src="PCA_analysis_files/figure-markdown_github/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

``` r
p3 <-  ggplot(dataGG, aes(PC1, PC2, color= condition, label = condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
p3 <- p3 + geom_text_repel(box.padding = 0.1) + scale_radius(range = c(3,6))
p3
```

<img src="PCA_analysis_files/figure-markdown_github/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

### Reference

Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8

Pham, V.N., Xu X, Huq E. Molecular bases for the constitutive photomorphogenic phenotypes in Arabidopsis. Development 2018:dev.169870 doi: 10.1242/dev.169870.
