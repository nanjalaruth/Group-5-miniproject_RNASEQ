#!/bin/bash
#PHASE1


#QUALITY CONTROL

#for sample in `cat sample_id.txt`
#do
#fastqc data/${sample}_R1.fastq.gz data/${sample}_R2.fastq.gz --outdir Fastqc_reports
#done
#multiqc Fastqc_reports -o  Fastqc_reports

#for sample in `cat sample_id.txt`
#do
      # trim_galore -j 10 -paired data/${sample}_R1.fastq.gz data/${sample}_R2.fastq.gz -q 25 --length 20 --fastqc -o Trim_galore
#done

#multiqc Trim_galore -o Trim_galore

#########################################################################################################################################

#PHASE2

#ALIGNMENT USING STAR

#CREATING A GENOME INDEX

#STAR --runThreadN 20 \
#--runMode genomeGenerate \
#--genomeDir /home/fkakembo/Bioinformatics/RNASeq_MiniProject/Ref_Genome \
#--genomeFastaFiles /home/fkakembo/Bioinformatics/RNASeq_MiniProject/References/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#--sjdbGTFfile /home/fkakembo/Bioinformatics/RNASeq_MiniProject/References/Homo_sapiens.GRCh38.100.gtf \
#--sjdbOverhang 99



#ALIGNING READS

#mkdir -p results/STAR
#for sample in `cat sample_id.txt`
#do

#STAR --genomeDir /home/fkakembo/Bioinformatics/RNASeq_MiniProject/Ref_Genome \
#--runThreadN 20 \
#--readFilesIn Trim_galore/${sample}_R1_val_1.fastqc Trim_galore/${sample}_R2_val_2.fastqc  \
#--outFileNamePrefix /results/STAR/${sample}. \
#--outSAMtype BAM SortedByCoordinate \
#--outSAMattributes Standard
#done

#########################################################################################################################################################


#ALIGNMENT USING HISAT2

##HISAT2

#hisat2-build -p 10 /home/fkakembo/Bioinformatics/RNASeq_MiniProject/References/Homo_sapiens.GRCh38.dna.primary_assembly.fa reference
#for sample in `cat sample_id.txt`
#do
#hisat2  -p 20 --summary-file ${sample}.summary.txt -x Hisat/reference -1 /home/sdlamini/Miniproject/Trim_galore/${sample}_R1_val_1.fq -2 /home/sdlamini/Miniproject/Trim_galore/${sample}_R2_val_2.fq -S ${sample}.sam
#done


#for sample in `cat sample_id.txt`
#do
#samtools view -Sb Sam/${sample}.sam | samtools sort  > ${sample}.sorted.bam
#samtools index ${sample}.sorted.bam
#done

###########################################################################################################################################################



#FEATURE COUNTS

#mkdir featurecounts
#for sample in `cat sample_id.txt`
#do
    #featureCounts -T 4   \
  # -a /home/fkakembo/Bioinformatics/RNASeq_MiniProject/References/Homo_sapiens.GRCh38.100.gtf  \
  # -o featurecounts/${sample}.txt \
  # Bam/${sample}.sorted.bam
#done


#GENERATING A MATRIX AND EXTRACTING READ COUNTS ASSOCIATED TO GENES

#mkdir Matrix
#for sample in `cat sample_id.txt`
#do
#cut -f1,7,8,9,10,11 featurecounts/${sample}.txt > Matrix/${sample}.Rmatrix.txt                                                                                         #done

#view this output
#less results/counts/sample37.Rmatrix.txt

#getting some stats of mapped and unmapped reads
#for sample in `cat sample_id.txt`
#do
 #  samtools flagstat Sam/${sample}.sam > ${sample}.txt
#done



#KALLISTO PSEUDO-ALIGNMENT


#building an index#
#kallisto index -i human /home/fkakembo/Bioinformatics/RNASeq_MiniProject/References/gencode.v34.transcripts.fa 
#for sample in `cat sample_id.txt`
#do 
   # kallisto quant -t 20 \
       # -i human -o KALI/${sample} \
    #Trim_galore/${sample}_R1_val_1_fastqc.zip \
    #Trim_galore/${sample}_R2_val_2_fastqc.zip
#done    


####################################################################################################################################################

#MULTIQC
#multiqc featurecounts
#multiqc summary
