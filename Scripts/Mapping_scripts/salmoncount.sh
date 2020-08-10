salmon index -t transcriptome/gencode.v34.transcripts.fa -i salmonindex/salmon_index

for sample in `cat sample_id.txt`

do 
        salmon quant -l A -i salmonindex/salmon_index -p 8 \
                -1 ../trimmed_reads/${sample}_R1_paired.fq \
                -2 ../trimmed_reads/${sample}_R2_paired.fq \
                --validateMappings -o salmon_output/${sample}_quant

done

