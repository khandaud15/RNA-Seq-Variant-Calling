shell.prefix("source ~/.bash_profile; set -euo pipefail;")
from util.varsub import varsub
configfile: "config.yaml"
varsub(config)
import glob
import os

SAMPLES, = glob_wildcards(config['datadirs']['fastq'] + "/" + "{file}_1.fq.gz")
READS = ["1", "2"]



rule all:
      input:
         config['reference']['stargenomedir']['hg38'] + "/" + "SAindex"
         config['datadirs']['sj_files'] + "/" + "SJ.out.pass1_merged.tab",
         expand(config['datadirs']['qc'] + "/" + "{file}_{read}_fastqc.html", file=SAMPLES, read= READS),
         expand(config['datadirs']['trim'] + "/" + "{file}_{read}_val_{read}.fq.gz", file = SAMPLES, read= READS),
         expand(config['datadirs']['bam'] + "/" + "{file}_SJ.out.tab", file = SAMPLES),
         expand(config['datadirs']['pass2'] + "/" + "{file}_Aligned.sortedByCoord.out.bam", file = SAMPLES ),
         expand(config['datadirs']['RGbam'] + "/" + "{file}_Aligned.sortedByCoord.out.RG.bam", file=SAMPLES),
         expand(config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam",file=SAMPLES),
         expand(config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam.bai", file=SAMPLES),
         expand(config['datadirs']['splitNcigar'] + "/" + "{file}_split.out.bam", file=SAMPLES),
         expand(config['datadirs']['Recal1'] + "/" + "{file}_recal.table", file=SAMPLES),
         expand(config['datadirs']['BQSR_1'] + "/" + "{file}_recal.pass1.bam", file=SAMPLES )


rule fastqc:
    input:
	     r1 = config['datadirs']['fastq'] + "/" + "{file}_1.fq.gz",
       r2 = config['datadirs']['fastq'] + "/" + "{file}_2.fq.gz"
    output:  config['datadirs']['qc'] + "/" + "{file}_{read}_fastqc.html"
    params:
	     prefix =  config['datadirs']['qc'], 
	  resources:
	     mem_mb= 10000
    threads : 12
    shell:"""
         fastqc  --threads {threads} --outdir {params.prefix} --nogroup {input.r1} {input.r2} 
         """

rule trim_galore_pe:
    input:
      f1 = config['datadirs']['fastq'] + "/" + "{file}_1.fq.gz",
      f2 = config['datadirs']['fastq'] + "/" + "{file}_1.fq.gz"
    output:
      fwd_pai = config['datadirs']['trim'] + "/" + "{file}_1_val_1.fq.gz",
      rev_pai = config['datadirs']['trim'] + "/" + "{file}_2_val_2.fq.gz",
    params:
      extra = " -j 8 --illumina -q 20 --phred33 --length 20",  
      prefix =  config['datadirs']['trim'],
    resources:
      mem_mb= 20000
    shell:""" 
        trim_galore \
        {params.extra} \
        --paired {input.f1} {input.f2} \
        -o {params.prefix} \
        --fastqc
         """

rule pass1:
   input:
      f1 = config['datadirs']['trim'] + "/" + "{file}_1_val_1.fq.gz",
      f2 = config['datadirs']['trim'] + "/" + "{file}_2_val_2.fq.gz",
      queue = rules.trim_galore_pe.output.rev_pai
   output: config['datadirs']['bam'] + "/" + "{file}_SJ.out.tab", config['datadirs']['bam'] + "/" + "{file}_Aligned.toTranscriptome.out.bam"
   params:
      genomedir = config['reference']['star_ref'],
      prefix =  config['datadirs']['bam'] + "/" + "{file}_"
   threads: 16
   resources:
      mem_mb= 40000
   shell: """
        STAR  \
        --runThreadN {threads} \
        --genomeDir {params.genomedir} \
        --readFilesIn {input.f1} {input.f2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype None \
        --outSAMunmapped Within \
        --quantMode TranscriptomeSAM \
        --outSAMattributes NH HI AS NM MD \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --limitBAMsortRAM 50000000000
        """ 


 
rule SJ_Merge:
    input:
      sjs =  expand(config['datadirs']['bam'] + "/" + "{file}_SJ.out.tab" , file = SAMPLES)
    output:
      sjs=  config['datadirs']['sj_files'] + "/" + "SJ.out.pass1_merged.tab"
    threads: 1
    shell: """
         cat {input.sjs} | awk '$7 >= 3' | cut -f1-4 | sort -u > {output.sjs}
          """

 
rule star_genome:
    input:
        fasta = config['reference']['fasta']['hg38'],
        gtf = config['reference']['gtf']['hg38'],
        sjs =  config['datadirs']['sj_files'] + "/" + "SJ.out.pass1_merged.tab",
        genomedir = config['reference']['stargenomedir']['hg38'],
        queue = rules.merge.output.sjs
    output:
        starindex = config['reference']['stargenomedir']['hg38'] + "/" + "SAindex"
    params:
        overhang = 99
    threads: 12
    resources:
        mem_mb = 40000
    shell: """
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {input.genomedir} \
        --outFileNamePrefix {input.genomedir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --limitSjdbInsertNsj 2037800 \
        --sjdbFileChrStartEnd  {input.sjs} \
        --sjdbOverhang {params.overhang}
          """

          
rule pass2:
   input:
      f1 = config['datadirs']['trim'] + "/" + "{file}_1_val_1.fq.gz",
      f2 = config['datadirs']['trim'] + "/" + "{file}_2_val_2.fq.gz",
      line = rules.star_genome.output.starindex
   output: config['datadirs']['pass2'] + "/" + "{file}_Aligned.toTranscriptome.out.bam", config['datadirs']['pass2'] + "/" + "{file}_Aligned.sortedByCoord.out.bam"
   params:
      genomedir = config['reference']['stargenomedir']['hg38'],
      prefix =  config['datadirs']['pass2'] + "/" + "{file}_"
   threads: 16
   resources:
      mem_mb= 50000
   shell: """
        STAR  \
        --runThreadN {threads} \
        --genomeDir {params.genomedir} \
        --readFilesIn {input.f1} {input.f2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --quantMode TranscriptomeSAM \
        --outSAMattributes NH HI AS NM MD \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --outBAMsortingThreadN 5 \
        --limitBAMsortRAM 50000000000
        """


rule AddRG:
     input:
         bam = config['datadirs']['pass2'] + "/" + "{file}_Aligned.sortedByCoord.out.bam"
     output:
         RG = config['datadirs']['RGbam'] + "/" + "{file}_Aligned.sortedByCoord.out.RG.bam"
     params: "RGLB=lib1 RGPL=illumina RGPU={file} RGSM={file}"   
     shell:"""
          module load picard
          java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups {params} I={input.bam} O={output.RG}
           """ 
                  

rule mark_dups:
    input:
       bam = config['datadirs']['RGbam'] + "/" + "{file}_Aligned.sortedByCoord.out.RG.bam"
    output:
       dbam = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam",
       metric = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.metrics.txt"
    params:
       picard = "java -jar $EBROOTPICARD/picard.jar"
    resources:
       mem_mb = 10000
    shell: """
         module load picard
        {params.picard} MarkDuplicates  INPUT={input.bam} OUTPUT={output.dbam} METRICS_FILE={output.metric} ASSUME_SORT_ORDER=coordinate OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
          """

 rule index:
      input:
         bam = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam"
      output:
         bai = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam.bai"
      shell:"""
            module load samtools
            samtools index {input.bam} {output.bai}
            """          

rule splitNcigar:
     input:
        bam = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam",
        fasta = config['reference']['fasta']['hg38']
     output:
        SBam = config['datadirs']['splitNcigar'] + "/" + "{file}_split.out.bam"
     resources:
        mem_mb= 60000
     shell:""" 
           module load gatk     
           gatk SplitNCigarReads \
           -R {input.fasta} \
           -I {input.bam} \
           -O {output.SBam} 
           """ 


rule BQSR_Pass1          
     input:
        bam = config['datadirs']['splitNcigar'] + "/" + "{file}_split.out.bam",
        1000G_SNPs = config['reference']['1000G']['hg38'],
        Indels = config['reference']['Indels']['hg38'],
        DbSNP = config['reference']['DbSNP']['hg38'],
        fasta = config['reference']['fasta']['hg38']
     output:
        Recall =  config['datadirs']['Recal1'] + "/" + "{file}_recal.table"
     resources:
        mem_mb = 50000
     shell:"""
           gatk BaseRecalibrator \
           -I {input.bam} \
           -R {input.fasta} \
           --known-sites  {input.1000G_SNPs} \
           --known-sites  {input.Indels}  \
           --known-sites  {input.DbSNP} \
           -O {output.recall}
           """

rule ApplyBQSR:
     input:
        bam = config['datadirs']['splitNcigar'] + "/" + "{file}_split.out.bam",
        fasta = config['reference']['fasta']['hg38'],
        recal = config['datadirs']['Recal1'] + "/" + "{file}_recal.table"
    output:
        Rbam = config['datadirs']['BQSR_1'] + "/" + "{file}_recal.pass1.bam"
    resources:
        mem_mb = 50000
    shell:"""
          gatk ApplyBQSR \
          -I {input.bam}  \
          -R {input.fasta} \
          --bqsr-recal-file {input.recal} \
          -O {output.Rbam}
          """