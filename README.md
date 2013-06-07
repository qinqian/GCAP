## What is GCAP
> `G`lobal `C`hromatin `A`ccessibility `P`ipeline
> Dnase, Bnase, Cnase的pipeline


Install 
=============

First, install pipeline framework: 

	hg clone https://bitbucket.org/hanfeisun/samflow
	# or
	pip install samflow
	
	
Install FastQC
------------------
*Site* 	<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
This part will be replaced by Jim's codes

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

No matter it's SE or PE data, we take 5' tags for hotspot peaks calling.

NOTE: 
	
	SPOT score is calculated by *-both-passes/*hotspot.twopass.zscore.wig
	*-both-passes/*.twopass.merge150.wgt10.zgt2.wig minimally thresholded hotspots
	*-both-passes/*.hotspot.twopass.fdr0.05.merge.[wig/bed] FDR thresholded hotspots

	take the encode 2 mouse data as an example.
	
	The following encode_mouse_treat_rep1 is the prefix.
	a.  This is top 10 of the loose hotspot regions encode_mouse_treat_rep1.hotspot.twopass.zscore.wig( totally 303941 regions ), which is used for SPOT score calculation.chr1    3322171 3322200 2.148710
	chr1    3346446 3346569 2.896090
	chr1    3359049 3359208 3.776190
	chr1    3360959 3361032 2.148710
	chr1    3406163 3406545 33.391200
	chr1    3406686 3406964 2.148710
	chr1    3407779 3407880 2.148710
	chr1    3412749 3412880 2.148710
	chr1    3413164 3413314 2.773560
	chr1    3415304 3415416 2.258360
	
	b. This is top 10 of `hotspot peaks` regions(Hotspot, broad peaks) encode_mouse_treat_rep1.twopass.merge150.wgt10.zgt2.wig
	so-called minimally thresholded hotspots ( totally 197552 regions )
	chr1    3322171 3322200 2.148710
	chr1    3346446 3346569 2.896090
	chr1    3359049 3359208 3.776190
	chr1    3360959 3361032 2.148710
	chr1    3406163 3406964 33.391200
	chr1    3407779 3407880 2.148710
	chr1    3412749 3412880 2.148710
	chr1    3413164 3413314 2.773560
	chr1    3415304 3415416 2.258360
	chr1    3427261 3427386 2.938350
	
	c. This  is the top ten of two passes hotspot regions with fdr 0.01 
	encode_mouse_treat_rep1.hotspot.twopass.fdr0.01.bed ( totally 162840 regions ), so called FDR thresholded hotspots.
	chr1    3406163 3406545 33.3912
	chr1    3445667 3445915 8.39721
	chr1    3467832 3468094 12.7712
	chr1    3504561 3504909 106.159
	chr1    3504916 3505327 192.588
	chr1    3505337 3505540 44.0136
	chr1    3505610 3505694 8.31948
	chr1    3541336 3541619 7.77236
	chr1    3541699 3541978 5.27296
	chr1    3542027 3542248 9.02206
	
	d. If wavelets peaks is so little that chromosome check is not True, above c will be copied to d.(Peaks, narrow peaks)
	This is the top 10 of merged wavelets peaks and hotspot regions with fdr 0.01
	encode_mouse_treat_rep1.hotspot.twopass.fdr0.01.merge.pks.bed ( totally 122532 regions, the most strict one, this is what I choose as final peaks number ), so called FDR thresholded peaks.
	chr1    3406300 3406450 .       36.406250
	chr1    3445740 3445890 .       8.953125
	chr1    3467860 3468010 .       14.734375
	chr1    3505040 3505190 .       194.953125
	chr1    3541480 3541630 .       8.593750
	chr1    3542020 3542170 .       9.250000
	chr1    3543180 3543330 .       25.515625
	chr1    3576620 3576770 .       16.453125
	chr1    3595820 3595970 .       47.406250
	chr1    3601420 3601570 .       9.750000

	

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
For single end data, we used MACS2 predictd to predict fragment size and calculate standard deviation by using MACS2 *predict_model.R. Fragment size evaluation for PE and SE will be replaced by MACS2 in the future.

All default, `2G` memory, `4cpu` will be used.

----

install UCSC component
-------------------------
*Site* <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/>
`BedClip` is used to remove outlier chromosome locations.
`bedGraphToBigWig` is used to convert hotspot reads density files to bigwiggle.


With the help of Jim, we would use `bigWigCorrelate`, it's built-in gcap/pipeline-scripts/bigWigCorrelate(added to $PATH, `bedToBigBed` is needed) for replicates consistency evaluation on union peaks regions.


Built-in modules
-------------------
This part is `built-in modules from cistrome-application`.
include BedIO, FeatIO, Func.
Export pipeline-scripts/conservation_average.py to $PATH, needs `bx-python`.

Three modes of sampling:

- We use built-in function to do raw reads sampling from PE and SE FASTQ.
- Python function to sample reads from PE and SE SAM files， including `mappable and unmappable reads`.
- picard sampling for PE and SE SAM files. (May use many threads and memory)

We decide to use built-in to sample 5M mappable and unmappable reads from PE and SE SAM files and both transfer to hospot as single end data.


DHS
-------
Merged DHS from ENCODE narrow peaks is used as reference union DHS regions.



Install latex and jinja2
-------------------------
Check whether `pdflatex` is executable or not.
`jinja2` is a template module for rendering latex document:
	
	pip install jinja2
	# for options
	pip install argparse

- - - -


Usage
=============

refer to static/GCAP_pe.conf for pair end data, static/GCAP_se.conf for single end data.

If input is single end data, use `,` to separate replicates files.
If input is pair end data, use `,` to separate pairs, `;` to separate replicates.

`Input Format`
Only support fastq files now, bam files and reads bed files support will be added later.

`Keep duplicate`
is an important option for peaks calling. You could customize it by python conf files. To keep duplicates tags,
just `keep_dup = T` in `[hotspot]`.

instructions on conf files:
	
	[Basis]
	treat = rep_1_pair1, rep_1_pair2; rep_2_pair1, rep_2_pair2
	seq_type = pe   ## pe for pair end data, se for single end data
	user = qinq
	id = testid
	species = hg19
	output = results_path
	read_length = 50
	
	[tool]
	mapping = bowtie ## or bwa, for reads mapping
	peak_calling = hotspot ## or macs2, finish
	
	[picard]
	markdup = path   ## MarkDuplicates.jar path
	sort = path      ## SortSam.jar, for converting and sorting
	threads = 4    
	sample = path	 ## DownsampleSam.jar path
	
	[contaminate]    ## for library contamination evaluation, use bowtie for fast assessment only	hg19 = /mnt/Storage/data/Bowtie/hg19
	mm9 = /mnt/Storage/data/Bowtie/mm9
	rn4 = /mnt/Storage/home/qinq/projects/chilin/rn4_index/rn4
	
	[fastqc]
	threads = 5
	
	[bowtie] ## if using bowtie
	max_aign = 1  ## keep unique mappable reads or not
	
	[lib] ## external data
	genome_index = /mnt/Storage/data/Bowtie/hg19                                       ## bowtie or bwa index 
	chrom_bed =  /mnt/Storage/home/qinq/lib/chr_limit_hg19.bed                         ## chromosome limitation BED file
	dhs = /mnt/Storage/data/DHS/DHS_hg19.bed                                           ## union DHS sites
	velcro = /mnt/Storage/home/qinq/lib/wgEncodeHg19ConsensusSignalArtifactRegions.bed ## black list
	phast = /mnt/Storage/data/sync_cistrome_lib/conservation/hg19/placentalMammals/    ## Phastcon score from UCSC
	tss = /mnt/Storage/home/qinq/lib/refgenes/hg19.refgene.tss                         ## +- 1kb tss extracted from UCSC refgene 
	
	[hotspot]
	chrom_info = chromosome_info   ## from hotspot website	
	mappable_region = path  ## this is downloaded from hotspot website, it's a necessary part for evaluating genomic promotor percentage
	keep_dup = T            ## duplicates or not
	
	## if your peak caller is macs2	, fill the following parameters
	[macs]
	species = hs
	keep_dup = all ## keep all duplicate tags
	shiftsize = 50 ## could be customized
	
	
for human should be the data contained in the static/hg19.refgene.tss, which is +- 1Kb from refseq tss, you could use your customized `promotor` regions


* dry run

only print command line: 

	GCAP run -c GCAP_se.conf

* real run

run with real data: 

	GCAP run -c GCAP_pe.conf

* resume

resume process when problems occurs:

	GCAP run -c GCAP_pe.conf --resume


* Step control

	- skip steps:

    		GCAP run -c GCAP_pe.conf --skip 9 --resume

	- from and to:
	
		   	GCAP run -c GCAP_pe.conf --from 1 --to 3 --resume
		   	
* clean up and purge
  to remove temporary files, only keep QC report, bigwiggle, peaks, hotspot and bam files, json and latex files:
  
  		GCAP clean -c GCAP_pe.conf
  		GCAP purge -c GCAP_pe.conf


----

Prototype  Features
======================

1. estimate per sequence quality and library contamination by using 100k sampled reads(mapping by bowtie and bwa(not added yet))
2. reads mapping and peaks calling on replicates all reads and 5M sampled reads(peak calling by hotspot and MACS2(added), sampling by Picard DownSampling by probability(5M/total_reads), which seems to be strange, sampling has been replaced by built-in python function, this is an option in conf file. set `picard sample path` would choose picard sampling, other situation would use built-in. Considering sampling from pair 1 for PE fastq or SE fastq files.
3. For pair end data, 5' tags from each pair will be treated as single end for hotspot v3.
4. peaks calling on combo of all reads and 5M sampled readss, Hotspot for 5M reads, Peaks for all reads, that is, use `b, d`.
5. estimate library complexity/redundancy by 5M reads(picard, Markduplicates)
6. estimates SPOT score for 5M reads (Hotspot, a), optional: MACS2 spot score
7. estimate replicates consistency by
	1. BigWiggle Correlation on 5M sampled data union hotspot(Hotspot, b, which is filtered to remove blacklist and outlier regions) by bigwiggle(bigWigCorrelate, merged by bedops -m)
    2. Overlap by hotspot(Hotspot, b) regions overlap from 5M reads(Intersection over Union regions, bedtools)
    
8. This has been removed, calculate hotspot(filtered Hotspot, b) promotor percentage and compare with genome promotor percentage for 5M reads Hotspot(b) regions(built-in script)
9. calculate Phastcon score of top 1000 non-promotor Hotspot(filtered Hotspot, b)regions in 100 bp width around summits for 5M reads(modified cistrome built-in module)
10. estimate (narrow peaks, d) overlap with ENCODE narrow peaks union DHS on replicates of 5M reads(bedtools)

Optionally
===========

Install MACS2 for optional peaks caller and SE fragment size and standard deviation estimating: 

	git clone https://github.com/taoliu/MACS/
	or 
	pip install MACS2

`Add MACS2 SPOT score calculation, much stringent than Hotspot spot score.`

Install bwa for optional reads mapping tool. 

![QC report](example.png)
