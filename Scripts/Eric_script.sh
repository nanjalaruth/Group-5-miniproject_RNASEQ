#mkdir -p FastQC_out

#for sample in `cat sampleIDs.txt`
#do
#	echo $sample
#	fastqc data/${sample}_R1.fastq.gz data/${sample}_R2.fastq.gz --outdir FastQC_out
#done

#multiqc FastqQC -o  FastqQC
#for sample in `cat sampleIDs.txt`
#do
#	echo $sample
#	trim_galore -j 10 --paired data/${sample}_R1.fastq.gz data/${sample}_R2.fastq.gz -q 25 --length 20 --fastqc -o trim_galore

#done 
#PE -threads 8 -phred33 data/${sample}_R1.fastq.gz data/${sample}_R2.fastq.gz -q 25 ILLUMINACLIP:TruSeq3-SE:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
#multiqc trim_galore -o trim_galore


STAR --runThreadN  --genomeDir /home/fkakembo/Bioinformatics/RNASeq_MiniProject/References/ --readFilesCommand zcat sample37_R1_val_1.fq.gz,sample38_R1_val_1.fq.gz,sample39_R1_val_1.fq.gz,sample40_R1_val_1.fq.gz,sample41_R1_val_1.fq.gz,sample42_R1_val_1.fq.gz sample37_R2_val_2.fq.gz,sample38_R2_val_2.fq.gz,sample39_R2_val_2.fq.gz,sample40_R2_val_2.fq.gz,sample41_R2_val_2.fq.gz,sample42_R2_val_2.fq.gz
