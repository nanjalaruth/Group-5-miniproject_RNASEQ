#!/bin/bash
echo -e "This is kallisto script to perform pseudo-alignment\n"
#--------------------"Lets start by installing kallisto"----------------------------------
#conda install -c bioconda kallisto

"We shall perform kallisto using two approaches"
# using the reference Ensembl human transcriptome V96 to build a kallisto index
# Using the prebuild kallisto index

#-------------"Now building a kallisto index from a reference transcriptome"--------------
#" Lets first donload the human transcriptome reference"
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz

#----------------"Now we build a kallisto index"-------------------------------------------
#kallisto index gencode.v34.transcripts.fa.gz -i kallisto_index -k 31
#----------------" Now let's run the quantification algorithm"---------------------------
for sample in `cat ../sample_id`
do
	kallisto quant --bias -t 25 \
		-i kallisto_index -o ${sample} \
		trimmed_data/${sample}_R1_val_1.fq trimmed_data/${sample}_R2_val_2.fq
done
kallisto quant -i kallisto_index \
	-o quant_2 --bias -t 20 \
	trimmed_data/sample37_R1_val_1.fq trimmed_data/sample37_R2_val_2.fq \
	trimmed_data/sample38_R1_val_1.fq trimmed_data/sample38_R2_val_2.fq \
	trimmed_data/sample39_R1_val_1.fq trimmed_data/sample39_R2_val_2.fq \
	trimmed_data/sample40_R1_val_1.fq trimmed_data/sample40_R2_val_2.fq \
	trimmed_data/sample41_R1_val_1.fq trimmed_data/sample41_R2_val_2.fq \
	trimmed_data/sample42_R1_val_1.fq trimmed_data/sample42_R2_val_2.fq

#--------------" Now running the quantification algorithm using prebuilt kallisto index"---------------------------

#"Lets first download the prebuilt kallisto index"
#wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz
# "Then we decompress the downloaded file"
#tar xzvf homo_sapiens.tar.gz
#" Now quantification of the transcripts using prebuilt index"
for sample in `cat ../sample_id`
do
	kallisto quant -i homo_sapiens/transcriptome.idx \
		-o prebuilt_quant_1/${sample} --bias -t 20 \
		trimmed_data/${sample}_R1_val_1.fq trimmed_data/${sample_R2_val_2.fq
done
# Even here, the for loop with kallisto did not work, any one who managed to run kallisto with a for loop?
#OR
kallisto quant -i homo_sapiens/transcriptome.idx \
	-o prebuilt_quant_2 --bias -t 20 \
	trimmed_data/sample37_R1_val_1.fq trimmed_data/sample37_R2_val_2.fq \
	trimmed_data/sample38_R1_val_1.fq trimmed_data/sample38_R2_val_2.fq \
	trimmed_data/sample39_R1_val_1.fq trimmed_data/sample39_R2_val_2.fq \
	trimmed_data/sample40_R1_val_1.fq trimmed_data/sample40_R2_val_2.fq \
	trimmed_data/sample41_R1_val_1.fq trimmed_data/sample41_R2_val_2.fq \
	trimmed_data/sample42_R1_val_1.fq trimmed_data/sample42_R2_val_2.fq
