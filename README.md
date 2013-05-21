## What is GCAP
> `G`lobal `C`hromatin `A`ccessibility `P`ipeline

Install 
=============

First, install pipeline framework: 

	hg clone https://bitbucket.org/hanfeisun/samflow
	# or
	pip install samflow
	
	
Install FastQC
------------------
*Site* 	<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>

Install Bowtie
----------------
*Site* <http://bowtie-bio.sourceforge.net/index.shtml>

Add binary to PATH and download `hg19, mm9, rn4` index for library contamination evaluation

install Hotspot
----------------
* Hotspot default setting is used
* Download hotspot v3 from <http://www.uwencode.org/proj/hotspot/>, get the hotspot5, wavelets, wavePeaks into `$PATH`, you may need to compile from hotspot-deploy directory

* Download `CHROM FILE` and `_MAPPABLE_FILE_` for your species and reads length
* Extract default SPOT output of Hotspot v3


NOTE: 
	
	SPOT score is calculated by *-both-passes/*hotspot.twopass.zscore.wig
	*-both-passes/*.twopass.merge150.wgt10.zgt2.wig minimally thresholded hotspots
	*-both-passes/*.hotspot.twopass.fdr0.05.merge.[wig/bed] FDR thresholded hotspots

	Use two pass hotspot and peaks merge files, that is `*-both-passes/*.hotspot.twopass.fdr0.01.merge.pks.bed` as peaks
	

install bedops and bedtools
-----------------------------
* use bedtools to get peaks overlap with union DHS region
* use bedops to merge replicates peaks regions
* use bedtools to get non overlap regions

*Site* <http://code.google.com/p/bedops/> <http://code.google.com/p/bedtools/>

*** 


install picard
--------------------
Use picard for SortSam, Markduplicates for both single end and pair end data.
Use picard for pair end data `median fragment size` and `fragment standard deviation` evaluation.
For single end data, we used MACS2 predictd

----

Built-in modules
-------------------
For this part, install cistrome-application temporarily. This will be replaced by `Bedtools` and `built-in modules from cistrome-application`.

Use `mercurial` to get necessary packages:

    hg clone https://bitbucket.org/cistrome/cistrome-applications-harvard
    
install each of the packages contained in `cistrome-applications-harvard` repo

* download gene_tables for CEAS and Phastcon score bigwiggles of Placetalmammals for human and mouse
* use the lastest repo


Install latex
---------------
Check whether `pdflatex` is executable or not.

- - - -


Usage
=============

refer to static/GCAP_pe.conf for pair end data, static/GCAP_se.conf for single end data.

If input is single end data, use `,` to separate replicates files.
If input is pair end data, use `,` to separate pairs, `;` to separate replicates.

instructions on conf files:
	
	[Basis]
	...
	treat = rep_1_pair1, rep_1_pair2; rep_2_pair1, rep_2_pair2
	...

	[lib] 
	tss = 
	...
	
for human should be the data contained in the static/hg19.refgene.tss, which is +- 1Kb from refseq tss, you could use your customized `promotor` regions


* dry run

only print command line: 

	GCAP.py run -c GCAP_se.conf

* real run

run with real data: 

	GCAP.py run -c GCAP_pe.conf

* resume

resume process when problems occurs:

	GCAP.py run -c GCAP_pe.conf --resume


* skip steps
GCAP.py run -c GCAP_pe.conf --skip 9 --resume

* from and to

GCAP.py run -c GCAP_pe.conf --from 1 --to 3 --resume

----

Prototype  Features
======================

1. estimate per sequence quality and library contamination by using 100k sampled reads

2. reads mapping and peaks calling on replicates all reads and 5M sampled reads

3. peaks calling on combo of all reads and 5M sampled reads

4. estimates SPOT score for 5M reads

5. estimate replicates consistency by using 5M sampled data union peaks by bigwiggle and by peaks regions overlap

6. calculate peaks promotor percentage and compare with genome promotor percentage for 5M reads

7. calculate Phastcon score of top 1000 non-promotor peaks regions in 100 bp width around summits for 5M reads

8. estimate peaks overlap with DHS on replicates and combo by 5M reads

Optionally
===========

Install MACS2 for optional peaks caller: 

	git clone https://github.com/taoliu/MACS/
	or 
	pip install MACS2
	
Install bwa for optional reads mapping tool.

![QC report](example.png)
