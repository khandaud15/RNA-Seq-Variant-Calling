# @author: Daud Khan
# @date: September 15, 2019
# @desc: Snakemake pipeline for https://www.broadinstitute.org/gatk/guide/article?id=3891

shell.prefix("source ~/.bash_profile; set -euo pipefail;")
from util.varsub import varsub
configfile: "config.yaml"
varsub(config)
import glob, os


# A snakemake regular expression matching the forward mate FASTQ files.
SAMPLES, = glob_wildcards(config['datadirs']['fastq'] + "/" + "{file}_1.fq.gz")

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
READS = ["1", "2"]

# List of "{sample}.g.vcf.gz"
# used for rule "combineGVCFs"

gvcfLst = expand(config['datadirs']['vcf'] + "/" + "{file}.g.vcf" , file=SAMPLES)



# Rules --------------------------------------------------------------------------------
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
         expand(config['datadirs']['BQSR_1'] + "/" + "{file}_recal.pass1.bam", file=SAMPLES),
         expand(config['datadirs']['Recal2'] + "/" + "{file}_recal.table", file=SAMPLES),
         expand(config['datadirs']['BQSR_2'] + "/" + "{file}_recal.pass2.bam", file=SAMPLES),
         expand(config['datadirs']['vcf'] + "/" + "{file}.g.vcf" , file=SAMPLES),
         config['datadirs']['CombinedGvcfs'] + "/" + "all.g.vcf"



# QC of raw fastq files.
rule fastqc:
   input:
      f1 = config['datadirs']['fastq'] + "/" + "{file}_{read}.fq.gz"
   output: config['datadirs']['qc'] + "/" + "{file}_{read}_fastqc.html", config['datadirs']['qc'] + "/" + "{file}_{read}_fastqc.zip"
   params:
      prefix =  config['datadirs']['qc'],
   resources:
      mem_mb= 10000
   shell:
        """
        fastqc  --thread 8 --outdir {params.prefix} --nogroup {input.f1}
        """

#Trimmming of the illumina adapters.
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
         
# 1. Map paired-end RNA-seq reads to the genome.
# 2. Count the number of reads supporting each splice junction.
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


# Merge the Splice junction informtaion from Pass1 Mapping  
rule SJ_Merge:
    input:
      sjs =  expand(config['datadirs']['bam'] + "/" + "{file}_SJ.out.tab" , file = SAMPLES)
    output:
      sjs=  config['datadirs']['sj_files'] + "/" + "SJ.out.pass1_merged.tab"
    threads: 1
    shell: """
         cat {input.sjs} | awk '$7 >= 3' | cut -f1-4 | sort -u > {output.sjs}
          """

# Make an index of the genome for STAR using the merged splice junction information to get better alignments around novel splice junctions.
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
        overhang = 149
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

# 1. Map paired-end RNA-seq reads to the genome.
# 2. Make a coordinate sorted BAM with genomic coordinates.
# 3. Count the number of reads mapped to each gene.
# 4. Count the number of reads supporting each splice junction.         
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

# add read groups,Platform
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
                  
# mark duplicates
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

# Index bam file using samtools
 rule index:
      input:
         bam = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam"
      output:
         bai = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam.bai"
      shell:"""
            module load samtools
            samtools index {input.bam} {output.bai}
            """ 


#Splits N Cigar Reads from bam file
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

# base recalibration
rule BQSR_Pass1:         
     input:
        bam = config['datadirs']['splitNcigar'] + "/" + "{file}_split.out.bam",
        GSNPs = config['reference']['1000G']['hg38'],
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
           --known-sites  {input.GSNPs} \
           --known-sites  {input.Indels}  \
           --known-sites  {input.DbSNP} \
           -O {output.Recall}
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

#Base Recalibration 
rule BQSR_Pass2:         
     input:
        bam = config['datadirs']['BQSR_1'] + "/" + "{file}_recal.pass1.bam",
        GSNPs = config['reference']['1000G']['hg38'],
        Indels = config['reference']['Indels']['hg38'],
        DbSNP = config['reference']['DbSNP']['hg38'],
        fasta = config['reference']['fasta']['hg38']
     output:
        Recall =  config['datadirs']['Recal2'] + "/" + "{file}_recal.table"
     resources:
        mem_mb = 50000
     shell:"""
           gatk BaseRecalibrator \
           -I {input.bam} \
           -R {input.fasta} \
           --known-sites  {input.GSNPs} \
           --known-sites  {input.Indels}  \
           --known-sites  {input.DbSNP} \
           -O {output.Recall}
           """ 

#detects systematic errors made by the sequencer when it estimates the quality score of each base call
rule ApplyBQSR:
     input:
        bam = config['datadirs']['BQSR_1'] + "/" + "{file}_recal.pass1.bam",
        fasta = config['reference']['fasta']['hg38'],
        recal = config['datadirs']['Recal2'] + "/" + "{file}_recal.table"
    output:
        Rbam = config['datadirs']['BQSR_2'] + "/" + "{file}_recal.pass2.bam"
    resources:
        mem_mb = 50000
    shell:"""
          gatk ApplyBQSR \
          -I {input.bam}  \
          -R {input.fasta} \
          --bqsr-recal-file {input.recal} \
          -O {output.Rbam}
          """  


#Variant Calling 
rule gatk_HaplotypeCaller:
    input:
        bam = config['datadirs']['BQSR_2'] + "/" + "{file}_recal.pass2.bam"
        fasta = config['reference']['fasta']['hg38']
    output:
        vcf =  config['datadirs']['vcf'] + "/" + "{file}.g.vcf"   
    resources:
        mem_mb = 50000
    shell:"""
           gatk HaplotypeCaller \
           -R {input.fasta} \
           -I {input.bam} \
           -ERC GVCF --output-mode EMIT_ALL_CONFIDENT_SITES  \
           --dont-use-soft-clipped-bases \
           -stand-call-conf 20.0  \
           -O {output.vcf}
           """
           
 #Combine all the gVCFS for joint calling 
rule CombineGvfs:
     input:
        vcfs =  gvcfLst,
        fasta = config['reference']['fasta']['hg38']
     output:
        combined = config['datadirs']['CombinedGvcfs'] + "/" + "all.g.vcf" 
     params:
        lst = " --variant " .join(gvcfLst)
     resources:
        mem_mb = 100000
     shell:"""
           module load gatk
           export JAVA_TOOL_OPTIONS=-Xmx140g
           gatk CombineGVCFs \
           -R {input.fasta} \
           -O {output.combined} \
           --variant {params.lst}
           """         
  
 

