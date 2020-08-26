#Global variables
SAMPLES, = glob_wildcards("data/{sample}_R1.fastq.gz")
ref_seq = "data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ref_annot= "References/Homo_sapiens.GRCh38.100.gtf"
ref_tpme= "References/gencode.v34.transcripts.fa"
r1="data/{sample}_R1.fastq.gz"
r2="data/{sample}_R2.fastq.gz"

import os

#Target rule
rule all:
        input:
	        "Results/Fastqc_out/multiqc_report.html", 
		"Results/trim_galore/multiqc_report.html",
		"staoutput/STAR_counts.txt",
		"quants/quant.sf"
#___________________FASTQC STEP_____________________#
rule fastqc_1:
        input:
                read1=r1,
                read2=r2

        output:
                "Results/Fastqc_out/{sample}_R1_fastqc.html",
                "Results/Fastqc_out/{sample}_R1.fastqc.zip",
                "Results/Fastqc_out/{sample}_R2_fastqc.html"
        log:
                "logs/Fastqc_out/{sample}.log"

        shell:
                "fastqc --outdir Results/Fastqc_out {input} {input.read1} {input.read2} -t 20"
                "2> {log}"
#________________MULTI_QC 1 STEP_____________________#
rule multiqc_1:
    """
    Aggregate all FastQC reports into a MultiQC report.
    """
        input:
                html=expand("Results/Fastqc_out/{sample}_R1_fastqc.html", sample=SAMPLES)
        params:
                dir="Results"

        output:
                html = "Results/multiqc_report.html",
#                stats = "Results/multiqc_general_stats.txt"
                                     
        shell:
                "multiqc -f {params.dir} -o {params.dir}"
#________________QUALITY TRIMMING____________________#
rule trim_galore:
        input:
                read1=r1,
                read2=r2,
                html = "Results/Fastqc_out/{sample}_R1.fastqc.html"


        output:
                "Results/trim_galore/{sample}_R1_val_1.fq.gz",
                "Results/trim_galore/{sample}_R2_val_2.fq.gz",
                 r2_html="Results/trim_galore/{sample}_R2_val_2_fastqc.html",
		 r1_html="Results/trim_galore/{sample}_R1_val_1_fastqc.html"
        shell:
                "trim_galore -q 25 -j 8 --length 20 --fastqc --paired {input.read1} {input.read2} -o Results/trim_galore"
#_________________MULTI_QC 2 STEP___________________#
rule multiqc2:
        input:
                trihtml=expand("Results/trim_galore/{sample}_R1_val_1_fastqc.html", sample=SAMPLES)
        output:
                "Results/trim_galore/multiqc_report.html"
        shell:
                "multiqc -f {input} -o Results/trim_galore"

#__________________STAR Alignment__________________#
#rule index_star:
#    """
#    Index genomic DNA sequences using STAR.
#    """
#    input: fa = rules.download_ensembl_dna.output.fa,
#           gtf = rules.download_ensembl_gtf.output.gtf
#    output: Genome = 'references/Genome'
#    threads: 8
#    run:
#        res_dir = os.path.dirname(output.Genome)
#        shell("""
#              STAR --runMode genomeGenerate        \
#                   --runThreadN {threads}          \
#                   --genomeFastaFiles {input.fa}   \
#                   --sjdbGTFfile {input.gtf}       \
#                   --genomeDir {res_dir}/
#              """)

rule align_star:
    """
    Align sequencing reads using STAR.
    """
    input: fq1 = "Results/trim_galore/{sample}_R1_val_1.fq.gz",
           fq2 = "Results/trim_galore/{sample}_R2_val_2.fq.gz",
           idx = "Ref_Genome/Genome"
    output: bam = 'staoutput/{sample}.star.bam',
            sj = 'staoutput/{sample}.star.sj'
    threads: 16
    run:
        res_dir = os.path.dirname(output.bam)
        idx_dir = os.path.dirname(input.idx)
        shell("""
              STAR --runThreadN {threads}                 \
                   --runMode alignReads                   \
                   --readFilesCommand pigz -d -c          \
                   --outSAMtype BAM Unsorted              \
                   --genomeDir {idx_dir}/                 \
                   --outFileNamePrefix {res_dir}/         \
                   --readFilesIn {input.fq1} {input.fq2}
              mv {res_dir}/Aligned.out.bam {output.bam}
              mv {res_dir}/SJ.out.tab {output.sj}   
              """)

rule sort_bam:
    """
    Sort bam by coordinates using sambamba.
    """
    input: bam = 'staoutput/{sample}.star.bam'
    output: bam = 'staoutput/{sample}.sorted.bam'
    params: mem = '35G'
    threads: 16
    run:
        tmp_dir = os.path.dirname(output.bam)
        shell('sambamba sort --tmpdir {tmp_dir} -t {threads} -m {params.mem} -o {output.bam} {input.bam}')
      
#_____________________FEATURECOUNTS______________________#
rule featurecounts:
        input:
                files=expand('staoutput/{sample}.sorted.bam', sample=SAMPLES),
                anno=ref_annot
        output:
                "staoutput/STAR_counts.txt"
        params:
                threads=16
        shell:
                "featureCounts -T {params.threads} -a {input.anno} -o {output} {input.files}"

#_________________________SALMON_____________________________#
rule Salmon_index:
        input:
                ref_tpme
        output:
                "ref_Transcriptome"
        shell:
                "salmon index -t {input} -i {output}"

rule salmon_quant:
        input:
                tri1=expand("Results/trim_galore/{sample}_R1_val_1.fq.gz", sample=SAMPLES),
                tri2=expand("Results/trim_galore/{sample}_R2_val_2.fq.gz", sample=SAMPLES),
                salidx="ref_Transcriptome"
        output:
                "quants/quant.sf"
        params:
        threads:16

        shell:
                "salmon quant -i {input.salidx} -l A -1 {input.tri1} -2 {input.tri2} -p 16 --validateMappings -o {output}"
#_______________________KALLISTO____________________________#
#rule Kallisto_index:
#        input:
#                ref_tpme
#        output:
#                "kallis/kallisto_index"
#        shell:
#                "kallisto index {input} -i {output}"

#rule kallisto_alignment:
#        input:
#                tr_1=expand("Results/trim_galore/{sample}_R1_val_1.fq.gz", sample=SAMPLES),
#                tr_2=expand("Results/trim_galore/{sample}_R2_val_2.fq.gz", sample=SAMPLES),
#                kalidx="kallis/kallisto_index"
#        output:
#                expand("kallisto/{sample}", sample=SAMPLES)
#        params:
#                threads=16
                
#        shell:
#                "kallisto quant -t {params.threads} -b 100 -i {input.kalidx} -o {output} {input.tr_1} {input.tr_2}"
