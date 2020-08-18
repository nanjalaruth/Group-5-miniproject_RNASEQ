mkdir -p FastqQC_Reports

for sample in `cat sample_id.txt`
do
#	echo $sample
	fastqc data/${sample}_R1.fastq.gz data/${sample}_R2.fastq.gz --outdir FastqQC_Reports
done

multiqc FastqQC_Reports -o  FastqQC_Reports

for sample in `cat sample_id.txt`
do
#	echo $sample
	trim_galore -j 10 --paired data/${sample}_R1.fastq.gz data/${sample}_R2.fastq.gz -q 25 --length 20 --fastqc -o Trim_galore

done 

multiqc Trim_galore -o Trim_galore

echo -e "\n Now Running Alignment \n"

#References/Homo_sapiens.GRCh38.100.gtf References/Homo_sapiens.GRCh38.dna.primary_assembly.fa  Ref_Genome/

STAR --runThreadN 15 \
	--runMode genomeGenerate \
	--genomeDir Ref_Genome \
	--genomeFastaFiles  References/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--sjdbGTFfile  References/Homo_sapiens.GRCh38.100.gtf 

for sample in `cat sample_id.txt`
do
	STAR  --runThreadN 25 \
		--genomeDir  Ref_Genome  \
		--readFilesIn  Trim_galore/${sample}_R1_val_1.fq.gz Trim_galore/${sample}_R2_val_2.fq.gz  \
		--readFilesCommand gunzip -c  \
		--outFileNamePrefix  STAR_Alignment/${sample}. 

#		--outSAMtype BAM SortedByCoordinate

	samtools view -Sb STAR_Alignment/${sample}.Aligned.out.sam  | samtools sort  > STAR_Alignment/${sample}_sorted.bam
	samtools index STAR_Alignment/${sample}_sorted.bam

	rm  STAR_Alignment/${sample}.Aligned.out.sam

done


##   USING HISAT2
hisat2-build -p 25 References/Homo_sapiens.GRCh38.dna.primary_assembly.fa   Hisat_Index/Homo_sapiens.GRCh38v3_hisat2.idx


for sample in `cat sample_id.txt`
do 
	hisat2 -p 20 -x Hisat_Index/Homo_sapiens.GRCh38v3_hisat2.idx  \
		-1 Trim_galore/${sample}_R1_val_1.fq.gz  -2 Trim_galore/${sample}_R2_val_2.fq.gz  \
		-S   HISAT_Alignment/${sample}_hisat.sam

	samtools view -Sb HISAT_Alignment/${sample}_hisat.sam  | samtools sort  > HISAT_Alignment/${sample}_hisat_sorted.bam
        samtools index HISAT_Alignment/${sample}_hisat_sorted.bam

	rm  HISAT_Alignment/${sample}_hisat.sam
done

##  Featurecounts

mkdir  -p  Featurecounts

#counts  for STAR
featureCounts -T 20 -a References/Homo_sapiens.GRCh38.100.gtf \
        -o Featurecounts/STAR_counts.txt \
	STAR_Alignment/sample{37..42}_sorted.bam


# counts for Hisat
featureCounts -T 20 -a References/Homo_sapiens.GRCh38.100.gtf \
	-o Featurecounts/hisat_counts.txt \
	HISAT_Alignment/sample{37..42}_hisat_sorted.bam



##  Pseudo alignment

# 1. Using kallisto
mkdir  -p  Kallisto

kallisto index References/gencode.v34.transcripts.fa -i Kallisto/kallisto_index

for sample in `cat sample_id.txt`
do
	kallisto quant -t 20 -b 100 \
		-i Kallisto/kallisto_index* -o Kallisto/${sample}   \
		Trim_galore/${sample}_R1_val_1.fq.gz Trim_galore/${sample}_R2_val_2.fq.gz
done


# 2. Using Salmon

salmon index -p 15 -i Salmon/salmon_index -t References/gencode.v34.transcripts.fa --gencode

for sample in `cat sample_id.txt`
do 
	salmon quant -l A -i Salmon/salmon_index -p 20  \
		-1 Trim_galore/${sample}_R1_val_1.fq.gz  \
		-2 Trim_galore/${sample}_R2_val_2.fq.gz  \
		--validateMappings -o Salmon/${sample}_quant
done
