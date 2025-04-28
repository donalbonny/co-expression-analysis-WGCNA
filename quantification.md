### Kallisto  https://pachterlab.github.io/kallisto/manual.html
# install kallisto
 https://pachterlab.github.io/kallisto/download
# download fasta file ensemble: https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index
kallisto quant -i Arabidopsis_TAIR10.idx  -o  output/${sample}  SRR3290788_S1_L001_R1_001.fastq.gz SRR3290788_S1_L001_R2_001.fastq.gz
# After quantification, kallisto generate the following file output 
abundance.h5
abudance.txt
run_info.json

# use h5 file and metadata file ( treatment) --> Heatmap and deseq2





