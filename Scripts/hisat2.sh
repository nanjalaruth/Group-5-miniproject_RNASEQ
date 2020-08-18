#!/bin/bash
echo -e "This is HISAT2  script"
# "Lets install hisat2 using conda"
#conda install -c bioconda hisat2
#for sample in `cat ../sample_id`
#do
#	hisat2 -p 25 -x Hisat_Index/Homo_sapiens.GRCh38v3_hisat2.idx   \
#		-1 ../trimmed_data/${sample}_R1_val_1.fq.gz -2 ../trimmed_data/${sample}_R2_val_2.fq.gz   \
#		-S Alignment_hisat/${sample}_hisat.sam

#	samtools view -Sb Alignment_hisat/${sample}_hisat.sam  | samtools sort  > Alignment_hisat/${sample}_hisat_sorted.bam
#	samtools index Alignment_hisat/${sample}_hisat_sorted.bam
done

echo -e " Using the FeatureCount to generate the transcript counts Hisat"
#featureCounts is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations
#It can be used to count both RNA-seq and genomic DNA-seq reads. It is available in the SourceForge Subread package or the Bioconductor Rsubread package.
# Here is the link "http://bioinf.wehi.edu.au/featureCounts/"

#To use featureCounts program included in the SourceForge Subread package, we have to first install subread package

#-------" Subread installation"----------------------------------------------------
# "conda installation"
# conda install -c bioconda subread

#---" Now generating the counts"
#mkdir -p FeatureCounts

#featureCounts -T 20 -a ../References/Homo_sapiens.GRCh38.100.gtf \ 
#	-o FeatureCounts/counts.txt \
#	Alignment_hisat/sample{37..42}_hisat_sorted.bam 

# tThis is the end of " Hisat2" script

