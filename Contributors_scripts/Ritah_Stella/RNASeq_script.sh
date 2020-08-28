
#Required in the working directory
# 1. sample_id.txt: contans a list of the sample names
# 2. Sample sequences (R1 and R2)

#==================================================================================================
#                                                                                                 #
#                       Data Preprocessing                                                        #
#                                                                                                 #
#==================================================================================================

echo -e "\n Data Preprocessing... \n"

mkdir -p QC_Reports  #create directory for the fastqc output

#Quality Check using 'FastQC'
for sample in `cat sample_id.txt` 
do
	fastqc ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz --outdir QC_Reports 
done

#Aggregate FastQC results using 'MultiQC'
multiqc QC_Reports -o  QC_Reports 

#Quality Trimming with 'TrimGalore'
mkdir -p Trim_galore

for sample in `cat sample_id.txt`
do
	trim_galore -j 10 --paired ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz -q 25 --length 20 --fastqc -o Trim_galore

# -j 10 specifies the number of cores to be used for trimming  
#--paired runs a paired-end validation on both trimmed R1 and R2 FastQ files to removes entire read pairs if at least one of the two sequences became shorter than the threshold
# -q 25 to remove portions of the reads with phred score less than 25 at the 3' end  
# --length 20 to filter trimmed reads less 20 bp

done 

#Quality Check the trimmed Reads with 'Mutliqc'
multiqc Trim_galore -o Trim_galore

#====================================================================================================
#                                                                                                   #
#                             Classical Alignment-based approach                                    #
#                                                                                                   #
#====================================================================================================

#Alignment of sequences to the reference genome using 'HISAT2'

echo -e "\n Now Running Classical Alignment-based with HISAT2... \n"

mkdir hisat2

# Download the human reference genome
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O hisat2/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip -d hisat2/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > hisat2/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Download the human annotaion file
wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz -p hisat2/Homo_sapiens.GRCh38.100.gtf.gz
gunzip -d hisat2/Homo_sapiens.GRCh38.100.gtf.gz > hisat2/Homo_sapiens.GRCh38.100.gtf

#HISAT2 indexing to build the reference genome index 
hisat2-build -p 25 hisat2/Homo_sapiens.GRCh38.dna.primary_assembly.fa hisat2/Homo_sapiens.GRCh38v3_hisat2.idx

#HISAT2 Alignment
for sample in `cat sample_id.txt`
do 
	hisat2 -p 20 -x hisat2/Homo_sapiens.GRCh38v3_hisat2.idx  \
		-1 Trim_galore/${sample}_R1_val_1.fq.gz  -2 Trim_galore/${sample}_R2_val_2.fq.gz  \
		-S   hisat2/${sample}_hisat2.sam

	samtools view -Sb hisat2/${sample}_hisat2.sam  | samtools sort  > hisat2/${sample}_hisat2_sorted.bam #compress the sam files to binary format
        samtools index hisat2/${sample}_hisat2_sorted.bam 

	rm  hisat2/${sample}_hisat2.sam #remove the .sam files
done

#Quantifying Aligned reads using the Subread package's script 'featureCounts'

featureCounts -T 20 -a hisat2/Homo_sapiens.GRCh38.100.gtf \
	-o hisat2/hisat2_counts.txt \
	hisat2/sample{37..42}_hisat2_sorted.bam

#====================================================================================================
#                                                                                                   #
#                            Kmer-based “Pseudo-alignment” approach                                 #
#                                                                                                   #
#====================================================================================================

# Aligning reads to the transcriptome using 'kallisto'

echo -e "\n Now Running Kmer-based “Pseudo-alignment” with Kallisto... \n"

mkdir -p Kallisto

# Download the human transcriptome reference
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz -P Kallisto/gencode.v34.transcripts.fa.gz

#Build the transcriptome index file
kallisto index Kallisto/gencode.v34.transcripts.fa.gz -i Kallisto/kallisto_index

for sample in `cat sample_id.txt`
do
	kallisto quant -t 20 -b 100 \
		-i Kallisto/kallisto_index* -o Kallisto/${sample}   \
		Trim_galore/${sample}_R1_val_1.fq.gz Trim_galore/${sample}_R2_val_2.fq.gz
done

echo -e "\n That's All folks!!! Go head with the Statistical Analyses in R \n"
