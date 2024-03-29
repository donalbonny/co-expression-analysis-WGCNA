# WGCNA analysis using RNA-seq datasets

My main codes follow the instruction from [Horvath lab](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html)

Useful links for the code 

1. Harvard DE Workshop  https://informatics.fas.harvard.edu/workshops/HarvardInformatics_DEworkshop_Fall2017.html

2. RNA-seq workshop at ASPB 2018, covered RNA-seq data analysis using Kallisto and Sleuth https://github.com/JasonJWilliamsNY/kallisto-rnaseq-jupyter/blob/master/sleuth.Rmd

3.  Kallisto https://pachterlab.github.io/kallisto/starting. Note that Kallisto is not used for mapping novel transcripts.
 
4. https://rafalab.github.io/dsbook/

5. UT Austin RNA-seq training  https://wikis.utexas.edu/display/bioiteam/Pseudomapping+with+Kallisto

6. I used MultiQC to examine the quality of multiple samples and give nice summary of FASTQC https://github.com/ewels/MultiQC 

7. UT Austin training on WCGNA  https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA

8. WGCNA slides https://pdfs.semanticscholar.org/ecac/5c1dc9dcbe6845abfd34e12ca091d3f993fb.pdf

9. Found the command line bootcamp which is really useful for beginner
http://rik.smith-unna.com/command_line_bootcamp/?id=tb3uo8zpxks

10. https://edu.isb-sib.ch/pluginfile.php/158/course/section/65/_01_SIB2016_wgcna.pdf

11. [Analysis of RNA-seq. PCA for outliers](https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html)


12. [Training on Differential gene expression, especially PCA to check sample quality](https://hbctraining.github.io/DGE_workshop/lessons/03_DGE_QC_analysis.html)

13. Lots of genomics data analysis references [here](http://www.begenomics.com/tutorial/#threeprnaseq) 

14. [co-expression tutorials and WGCNA understanding](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/)

15. [Bioconductor-RNA_seq analysis](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)

Get  gene name from biomaRt, Arabidopsis gene names. I figured out this is so much easier than the traditionally method that I used before to combine gene name from TAIR website with my results. 
 
Instead using biomaRt and choose athaliana_eg_gene  from host =plants.ensembl.org


 
 
Other turotials
* [Venn diagram](http://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html)

 

