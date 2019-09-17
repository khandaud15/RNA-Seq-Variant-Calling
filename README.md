## RNA-Seq Variant Calling Pipeline 
 This workflow is based on calling variants on RNA-Seq data using GATK4. the pipeline starts all the  way from raw Fastq files and end up with giving the VCF files. 

 ## Main Steps 

 #### Mapping to the Reference
 #### Tools involved:[STAR](https://github.com/alexdobin/STAR)
The pipeline begin with mapping RNA reads to a reference, we have used STAR aligner because it increased sensitivity compared to other alligner(especially for INDELS), as well as use STAR’s two-pass mode to get better alignments around novel splice junctions.

#### Add read groups, sort, mark duplicates, and create index using [Picard](https://broadinstitute.github.io/picard/) and [Samtools](http://www.htslib.org/doc/samtools.html)

The Star Mapping  step produces a BAM/SAM file, which we then put through the usual Picard processing steps: adding read group information, sorting, marking duplicates and indexing for downstream processing.

 #### Split'N'Trim and Reassign mapping qualities

This step  splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions as well reassign mapping qualities to the alligned reads because STAR Napping assigns good alignments a MAPQ of 255 (which technically means “unknown” and is therefore meaningless to GATK)