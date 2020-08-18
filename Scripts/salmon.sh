#!/bin/bash
echo -e "This is a salmon script for seudoalignment of RNA data\n"
# it involve creating and environment, installation,  building an index on a  transcriptome, and then quantifying RNASeq data
#conda create --name salmon #creating a salmon environment
#conda activate salmon
#wget https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz
#tar xzvf salmon-1.3.0_linux_x86_64.tar.gz

# This is the path to salmon after decompressing with the above command
#    "salmon-latest_linux_x86_64/bin/salmon"

echo -e "Obtaining Transcriptome and building an index\n" 

#-----------------"Now downloading the human transcriptome"---------------------------------------------
###### Here is the link to the reference transcriptome "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz"------------------------------------------------
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz

#-----------------"Now build an index on our transcriptome"---------------------------------------------
#salmon-latest_linux_x86_64/bin/salmon index -t gencode.v34.transcripts.fa.gz -i transcripts_index

#----------------"Now lets quantify the RNASeq data/samples"-------------------------------------------

for sample in `cat ../sample_id`
do
echo "Processing sample ${sample}"
	salmon-latest_linux_x86_64/bin/salmon quant -i transcripts_index -l A -1 ../trimmed_data/${sample}_R1_val_1.fq.gz -2 ../trimmed_data/${sample}_R2_val_2.fq.gz -p 20 --validateMappings -o quants/${sample}_quant
done
	
