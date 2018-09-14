# WGCNA analysis using RNA-seq datasets



Useful links for the code 

1. Harvard DE Workshop  https://informatics.fas.harvard.edu/workshops/HarvardInformatics_DEworkshop_Fall2017.html
2. RNA-seq workshop at ASPB 2018, covered RNA-seq data analysis using Kallisto and Sleuth https://github.com/JasonJWilliamsNY/kallisto-rnaseq-jupyter/blob/master/sleuth.Rmd
3.  Kallisto https://pachterlab.github.io/kallisto/starting
4. https://rafalab.github.io/dsbook/
5. UT Austin RNA-seq training  https://wikis.utexas.edu/display/bioiteam/Pseudomapping+with+Kallisto
6. I used MultiQC to examine the quality of multiple samples and give nice summary of FASTQC https://github.com/ewels/MultiQC 
7. UT Austin training on WCGNA  https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA


8. Found the command line bootcamp which is really useful for beginner
http://rik.smith-unna.com/command_line_bootcamp/?id=tb3uo8zpxks



# Get  gene name from biomaRt, Arabidopsis gene names. 

marts <- listMarts()
marts
marts <- listMarts(host = "plants.ensembl.org")
marts 
plants_mart <- useMart("plants_mart", host = "plants.ensembl.org" )
listDatasets(plants_mart)


plants_mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org" )
listAttributes(plants_mart)



t2g <- getBM(attributes = c("ensembl_transcript_id", 
                            "ensembl_gene_id", 
                            "description",
                            "external_gene_name"), 
             mart = plants_mart)
ttg <- dplyr::rename(t2g, target_id= ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
