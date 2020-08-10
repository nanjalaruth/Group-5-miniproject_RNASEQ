hisat2-build -p 8 ./reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./index_files/Homo_sapiens.idx

for sample in `cat sample_id.txt`

do

	hisat2 -p 8 -x ./index_files/Homo_sapiens.idx -1 ../trimmed_reads/${sample}_R1_paired.fq -2 ../trimmed_reads/${sample}_R2_paired.fq -S ./SAMfiles/${sample}.sam

	samtools view -b ./SAMfiles/${sample}.sam | samtools sort -o ./SortedBAM/${sample}_sorted.bam  

	samtools index ./SortedBAM/${sample}_sorted.bam

done 

