======
 Home
======

Section A: Sequence Quality
===========================

We use uniform FastQC_ program. We extracted the Phred score for
evaluation.

Section B: Library contamination
================================
We use bowtie to map each fastq file separately back to 3 species
genomes, mouse, rat and human. The right species should get the
highest mapping ratio.

Section C: Read mapping
=======================
bowtie for pair end mapping, use max insert size as 600 to get maximum
alignment. we trim the 3' with adapter size or default length.

Section D: Library complexity/redundancy
========================================
We need to convert sam files to bam in order to compress the file sizes.
samtools view -b
picard SamFormatConverter.jar

Before evaluating by picard, Bams should be sorted.
For pair end sequencing, We use picard tools for library complexity evaluation
::
   SortSam.jar I=testid.sam O=testid.sam.sorted SO=coordinate
   SortSam  26402145 10.56 minutes

picard uses so many cpus, default we could set to 4.
::
   MarkDuplicates.jar I=testid.sam.sorted
   O=testid_rep1.sam.sorted.markdup
   METRICS_FILE=testid_rep1.sam.metrics

For single end sequencing, we use all previous results as reference to
compare the newcoming data by house-made script.

Section E: Peak calling
=======================
Hotspot v3
MACS2 is an alternatively tool for DNase I peaks calling

Section F: Fragment sizes
=========================
picard:


macs2 predictd

Section G: Genome distribution(esp promoter %)
==============================================


Section H: Signal to noise
==========================
We use SPOT score to evaluate signal to noise ratio

* it is a by-product generated in hotspot pipeline, we just extract it
  and write into json file.

Section I: Overlap with historical DNase-seq
============================================


Section J: Evolutionary conservation
====================================


Section K: Replicate consistency
================================


* QC report summary information ::
     Give an overview of all the measurement pass or fail information

.. _FastQC site: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
