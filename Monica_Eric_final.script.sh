#!/bin/bash
echo -e "This script was written by Monica Mbabazi and Eric Kariuki\n"

#======================"RNASeq analysis workflow Bash script"===========================================================================
#This script covers the general workflow of RNASeq analysis for differential expression, and below are the steps:

# 1. "quality control and assurance of raw reads"
# 2. "Determining how many read counts are associated with known genes"
# we used two approaches  to determine number of read counts associated with known genes
#--------------1.'Classical Alignment-based approach'; aligning reads to a reference genome using 'Hisat2'
#--------------2.'Kmer-based “Pseudo-alignment” approaches'; for aligning reads to the transcriptome using 'kallisto'.
# 3. "Abandance Estimation " 
## This is only done for classical alignment based approach not applicable for pseudo-alignment, and used 'FeatureCounts' method
## After this step we had to import our "featurecounts" results and  "output" from kallisto into R to generate 'a DeSeq2' and 'EdgeR' objects for the downstream analysis.
#===========================================================================================================================
mkdir -p Results/fastqc

for sample in `cat sample_id`
do
echo -e "\n Performing Fastqc and quality control on $sample"
	fastqc data/${sample}* -o Results/fastqc

echo -e "\n performing quality control using Trim_galore\n"
	mkdir -p Results/trim_galore
	trim_galore -q 25 -j 20 --length 50 --fastqc -O Results/trim_galore --paired data/${sample}_R1.fastq.gz  data/${sample}_R2.fastq.gz 

done

echo -e "\n Running multiqc on fastqc results\n"
multiqc Results/fastqc -o Results/

echo -e "\n Running multiqc on trimmed results\n"
multiqc Results/trim_galore -o Results/trim_galore/

#=========================================================================================================================

echo -e "\n Determining how many read counts are associated with known genes\n"

#---------' Aligning reads to the reference genome'

# We used two methods to align the reads to the reference; 'STAR' and 'Hisat2' 

#========================================================================================================================
echo -e "\n Now running STAR"

#============= 'Installing STAR in my conda environment'
conda install -c bioconda star

#============= 'Aligning reads to the human genome reference using STAR'

#-----------------------------------------------------------------------------------------------------------------------------------
STAR --genomeDir /home/fkakembo/Bioinformatics/RNASeq_MiniProject/Ref_Genome --readFilesCommand zcat --readFilesIn sample37_R1_val_1.fq.gz,sample38_R1_val_1.fq.gz,sample39_R1_val_1.fq.gz,sample40_R1_val_1.fq.gz,sample41_R1_val_1.fq.gz,sample42_R1_val_1.fq.gz sample37_R2_val_2.fq.gz,sample38_R2_val_2.fq.gz,sample39_R2_val_2.fq.gz,sample40_R2_val_2.fq.gz,sample41_R2_val_2.fq.gz,sample42_R2_val_2.fq.gz
#============='"Running samtools on the .sam file"
for sample in `cat sample_id`
do
	samtools view -Sb staoutput/${sample}.sam | samtools sort  > ${sample}.sorted.bam
	samtools index ${sample}.sorted.bam
done
##-----------------------------------------------------------------------------------------------------------------------

# OR, also the for loop below can be used:
for sample in `cat sample_id`
do
	STAR  --runThreadN 25 \
		--genomeDir /home/fkakembo/Bioinformatics/RNASeq_MiniProject/Ref_Genome  \
		--readFilesIn  ${sample}_R1_val_1.fq.gz ${sample}_R2_val_2.fq.gz  \
		--readFilesCommand gunzip -c  \
		--outFileNamePrefix  staoutput/${sample}.

		--outSAMtype BAM SortedByCoordinate
	samtools view -Sb staoutput/${sample}.Aligned.out.sam  | samtools sort  > staoutput/${sample}_sorted.bam
	samtools index staoutput/${sample}_sorted.bam
done
#------------------------------------------------------------------------------------------------------------------------
echo "Multiqc of the STAR output"
multiqc staoutput/Log.final.out 
#========================================================================================================================

echo -e "\n Now running Hisat2"

#---------------------------- "Lets install hisat2 using conda in our conda environment"
# for this you need a human reference genome and then its index, or you can index the reference to generate the indexed reference
mkdir -p hisat2
cd hisat2
conda install -c bioconda hisat2
#-----'Building a reference genome index for Hisat2 usage'
hisat2-build -p 25 References/Homo_sapiens.GRCh38.dna.primary_assembly.fa   Hisat_Index/Homo_sapiens.GRCh38v3_hisat2.idx 

#----------' Running Hisat2 using the indexed reference'

for sample in `cat ../sample_id`
do
	hisat2 -p 25 -x Hisat_Index/Homo_sapiens.GRCh38v3_hisat2.idx   \
		-1 ../trimmed_data/${sample}_R1_val_1.fq.gz -2 ../trimmed_data/${sample}_R2_val_2.fq.gz   \
		-S Alignment_hisat/${sample}_hisat.sam

	samtools view -Sb Alignment_hisat/${sample}_hisat.sam  | samtools sort  > Alignment_hisat/${sample}_hisat_sorted.bam
	samtools index Alignment_hisat/${sample}_hisat_sorted.bam

	 rm  HISAT_Alignment/${sample}_hisat.sam
done
#==========================================================================================================================

#########################################################################################################################

echo -e 'Now using Kmer-based “Pseudo-alignment” approaches'

#---------"Aligning the reads to the reference transcriptome"------------------------------------------------------------
## Here we used to methods; 'Kallisto' and 'salmon' to do pseudo-alignment

#-----------------------------------------------------------------------------------------------------------------------
echo -e '\n Running salmon to align reads to the transcriptome'

#'Salmon' is a free software tool for estimating transcript-level abundance from RNA-seq read data.
#This involved salmon installation,  building an index on a  transcriptome, and then quantifying RNASeq data

#'Instead of installation, we obtained salmon from the pre-compiled banaries' using the code below;
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz

#Then decompress the downloaded file
tar xzvf salmon-1.3.0_linux_x86_64.tar.gz
# This is the path to salmon after decompressing with the above command
#    "salmon-latest_linux_x86_64/bin/salmon"
#-----------------"Now downloading the human transcriptome"---------------------------------------------
## Here is the link to the reference transcriptome "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz"
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz

#-----------------"Now build an index on our transcriptome"---------------------------------------------
salmon-latest_linux_x86_64/bin/salmon index -t gencode.v34.transcripts.fa.gz -i transcripts_index

#----------------"Now lets quantify the RNASeq data/samples"-------------------------------------------

for sample in `cat ../sample_id`
do
echo "Processing sample ${sample}"
        salmon-latest_linux_x86_64/bin/salmon quant -i transcripts_index -l A -1 ../trimmed_data/${sample}_R1_val_1.fq.gz -2 ../trimmed_data/${sample}_R2_val_2.fq.gz -p 20 --validateMappings -o quants/${sample}_quant
done
# Salmon outputs 6 files per sample but the most important file that contains transcript counts is 'quant.sf'
#---------------------------------------------------------------------------------------------------------------------------

echo -e '\n Running kallisto to align reads to the transcriptome'
#--------------------"Lets start by installing kallisto"----------------------------------
conda install -c bioconda kallisto

"We performed kallisto using two approaches"
#1. Using the reference Ensembl human transcriptome V96 to build a kallisto index
#2. Using the prebuild kallisto index

#-------------"Building a kallisto index from a reference transcriptome"--------------

#" Lets first donload the human transcriptome reference"
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz

#----------------"Now we build a kallisto index"-------------------------------------------
kallisto index gencode.v34.transcripts.fa.gz -i kallisto_index -k 31

#----------------" Now let's run the quantification algorithm"---------------------------
mkdir  -p  Kallisto

for sample in `cat ../sample_id`
do
        kallisto quant --bias -t 25 \
                -i kallisto_index -o ${sample} \
                trimmed_data/${sample}_R1_val_1.fq trimmed_data/${sample}_R2_val_2.fq
done
#------------------------------------------------------------------------------------------------------------------------

echo -e "\n Now running the quantification algorithm using prebuilt kallisto index"

#"Lets first download the prebuilt kallisto index"
wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz

# "Then we decompress the downloaded file"
tar xzvf homo_sapiens.tar.gz

#-------------------------" Now quantification of the transcripts using prebuilt index"--------------------------------
for sample in `cat ../sample_id`
do
        kallisto quant -i homo_sapiens/transcriptome.idx \
                -o prebuilt_quant_1/${sample} --bias -t 20 \
                trimmed_data/${sample}_R1_val_1.fq trimmed_data/${sample_R2_val_2.fq
done
##########################################################################################################################


#========================================================================================================================
#===================="GENERATING THE TRANSCRIPT COUNTS FROM HISAT2"======================================================
# This can be done using different methods like 'featureCount' and 'HTSeq'

echo -e " Using the featureCounts to generate the transcript counts Hisat"
#featureCounts is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations
#It can be used to count both RNA-seq and genomic DNA-seq reads. It is available in the SourceForge Subread package or the Bioconductor Rsubread package.
# Here is the link "http://bioinf.wehi.edu.au/featureCounts/"
#To use featureCounts program included in the SourceForge Subread package, we had to first install subread package

#-------------------------------------------" Subread installation"----------------------------------------------------
# "conda installation"

conda install -c bioconda subread

#-----------" Now generating the Hisat2 counts using faetureCounts method"

mkdir -p FeatureCounts

featureCounts -T 20 -a ../References/Homo_sapiens.GRCh38.100.gtf \ 
       -o FeatureCounts/counts.txt \
       Alignment_hisat/sample{37..42}_hisat_sorted.bam 

#-----------------------------------------------------------------------------------------------------------------------
echo -e " Using the featureCounts to generate the transcript counts from STAR output"
featureCounts -T 20 -a ../References/Homo_sapiens.GRCh38.100.gtf \
	-o Featurecounts/STAR_counts.txt \
	STAR_Alignment/sample{37..42}_sorted.bam
#==============================================================================================================================

# Now after generating the counts using the method of choice, you just import your counts and the metadata into R for statistical analysis.
#These counts can be used to do differential expression using different R packages like 'DeSeq2', 'EdgeR', and 'sleuth'.
################"This is end of our  bash script for RNASeq analysis"########################################################
