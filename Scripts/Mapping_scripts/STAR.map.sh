ls
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir RefGenome \
--genomeFastaFiles reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile reference_genome/Homo_sapiens.GRCh38.100.gtf
--sjdbOverhang 99

for sample in `cat sample_id.txt`
do
	STAR --runThreadN 8 \
	--genomeDir Ref_Genome \
	--readFilesIn  ../trimmed_reads/${sample}_R1_paired.fq ../trimmed_reads/${sample}_R2_paired.fq  \
	--outFileNamePrefix  ./STAR_Alignment/${sample}.

	samtools view -Sb ./STAR_Alignment/${sample}.Aligned.out.sam  | samtools sort  > ./STAR_Alignment/${sample}_sorted.bam

	samtools index ./STAR_Alignment/${sample}_sorted.bam
done
