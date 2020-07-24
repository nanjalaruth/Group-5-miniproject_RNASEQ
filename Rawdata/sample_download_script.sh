#! /bin/bash

for sample in `cat sample_id.txt`
do
	wget http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/${sample}_R1.fastq.gz
	wget http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/${sample}_R2.fastq.gz
done

