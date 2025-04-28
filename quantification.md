##### Kallisto  https://pachterlab.github.io/kallisto/manual.html
### install kallisto
 https://pachterlab.github.io/kallisto/download
### run quantification
kallisto quant -i ../../ref/Arabidopsis_TAIR10.idx  -o ../output/SRR1560618 SRR1560618_S1_L001_R1_001.fastq.gz  SRR1560618_S1_L001_R2_001.fastq.gz

### After quantification, kallisto generate the following file output 
abundance.h5
abudance.txt
run_info.json

### use h5 file and metadata file ( treatment) --> Heatmap and deseq2





