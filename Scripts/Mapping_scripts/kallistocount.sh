kallisto index -i kallindex/kallisto_index transcriptome/gencode.v34.transcripts.fa 

for sample in `cat sample_id.txt`

do

kallisto quant -t 8 -i kallindex/kallisto_index -o kallisto_output/${sample} ../trimmed_reads/${sample}_R1_paired.fq
 ../trimmed_reads/${sample}_R2_paired.fq

done
