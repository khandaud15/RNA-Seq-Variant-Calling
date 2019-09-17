## RNA-Seq Variant Calling Pipeline 
 This workflow is based on calling variants on RNA-Seq data using GATK4. the pipeline starts all the  way from raw Fastq files and end up with giving the VCF files. 

 ## Main Steps 

 #### Mapping to the Reference
 ###### Tools involved:[STAR](https://github.com/alexdobin/STAR)
The pipeline begin with mapping RNA reads to a reference, we have used STAR aligner because it increased sensitivity compared to other alligner(especially for INDELS), as well as use STAR’s two-pass mode to get better alignments around novel splice junctions.
