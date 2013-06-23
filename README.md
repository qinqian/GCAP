## What is GCAP
> `G`lobal `C`hromatin `A`ccessibility `P`ipeline
> 
> X-nase(Dnase, Bnase, Cnase) Quality analysis pipeline



### Installation 


NOTE:
	
	GCAP and samflow use python3
	
First, install pipeline framework: 

	
	hg clone https://bitbucket.org/hanfeisun/samflow
	python3 setup.py install
	# or
	pip-3.2 install samflow

Then, install GCAP:

    git clone https://github.com/qinqian/GCAP
    python3 setup.py install

The following component tools are all newest version.
	
##### Install FastQC
*Site* 	<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>, to keep `fastqc` in $PATH.
Info: 
	
	Per sequence quality score median is parsed from fastq text output by GCAP.

##### Install Bowtie and BWA

*Site* <http://bowtie-bio.sourceforge.net/index.shtml>

Add binary to PATH and download `hg19, mm9, rn4` index for bowtie and bwa mapping and library contamination evaluation.

Install bwa reads mapping tool, bwa is used for reads mapping evaluation only, which needs corresponding species index, default: 4 threads.

##### install Hotspot

* Hotspot default setting is used
* Download hotspot v3 from <http://www.uwencode.org/proj/hotspot/>, get the hotspot5, wavelets, wavePeaks into `$PATH`, you may need to compile from hotspot-deploy directory
* Download `CHROM FILE` and `_MAPPABLE_FILE_` for your species and reads length

No matter it's SE or PE data, we take 5' tags for hotspot peaks calling.

NOTE on hotspot output: 
	
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

	

##### install bedops and bedtools
* use bedtools to get peaks overlap with union DHS region
* use bedops to merge replicates peaks regions
* use bedtools to get non overlap regions(or bedops to filter blacklist)

*Site* <http://code.google.com/p/bedops/> <http://code.google.com/p/bedtools/>

*** 


##### install picard and samtools
- Use picard for SortSam, Markduplicates for both single end and pair end data, this will be replaced with `census`
- Use picard for pair end data `median fragment size` and `fragment standard deviation` evaluation.
- For single end data, we used MACS2 predictd to predict fragment size and calculate standard deviation by using MACS2 
  *predict_model.R.(this needs to be improved)

`2G` memory, `4cpu` will be used in picard Markduplicates and Insert size evaluation, `5G` memory and `4cpu` will be used for picard SortSam. If you want to use picard to sample large files to 5M, use -XX:-UseGCOverheadLimit and 4g memory setting.

----

##### install UCSC component
*Site* <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/>
`BedClip` is used to remove outlier chromosome locations.
`bedGraphToBigWig` is used to convert hotspot reads density files to bigwiggle.

Replicates consistency
a. With the help of Jim, we would use `bigWigCorrelate`, it's built-in gcap/pipeline-scripts/bigWigCorrelate(added to $PATH, `bedToBigBed` is needed) for replicates consistency evaluation on union DHS regions(filted by blacklist).
b. For whole genome correlation, use `wigCorrelate`


#### Built-in modules
Export gcap/pipeline-scripts/conservation_average.py to $PATH, needs `bx-python <https://bitbucket.org/james_taylor/bx-python/wiki/Home>` install.


#### DHS
Merged DHS from ENCODE narrow peaks(only wgEncodeUwDnase* narrowPeaks because of their high quality) is used as reference union DHS regions. We use ENCODE narrowPeak for union DHS extraction, we get union DHS from +/- 150bp from narrowPeak summits, *site* <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/> please email authors to get the BED files.

#### Blacklist
*Site* <http://compbio.tongji.edu.cn/~qinq/wgEncodeDacMapabilityConsensusExcludable.bed> is the latest blacklist for human, no mouse blacklist is obtained now.

##### fragment size tools

Install MACS2 for optional peaks caller and SE fragment size and standard deviation estimating: 

	git clone https://github.com/taoliu/MACS/
	or 
	pip install MACS2

`Add MACS2 SPOT score calculation, much stringent than Hotspot spot score.`


#### Install latex and jinja2
Check whether `pdflatex`(pdflatex (Version 3.141592-1.21a-2.2 (Web2C 7.5.4))) is executable or not. For d
`jinja2` is a template module for rendering latex document:
	
	pip-3.2 install jinja2
	# for options
	pip-3.2 install argparse

- - - -


### Usage

Three steps to use GCAP

##### 1. setup a python conf files

refer to gcap/static/GCAP_pe.conf for pair end data, gcap/static/GCAP_se.conf for single end data.
Example:
	
	[Basis]
	user = qinq
	id = testid
	species = hg19
	treat = rep_1_pair1, rep_1_pair2; rep_2_pair1, rep_2_pair2
	output = results_path
	read_length = 50
	sequence_type = pe
	
	[tool]
	mapping = bowtie ## or bwa, for reads mapping
	peak_calling = hotspot ## or macs2
	
	[picard]
	markdup = /mnt/Storage/home/qinq/softwares/picard-tools-1.91/MarkDuplicates.jar
	sort = /mnt/Storage/home/qinq/softwares/picard-tools-1.91/SortSam.jar
	threads = 4
	insertsize = /mnt/Storage/home/qinq/softwares/picard-tools-1.91/CollectInsertSizeMetrics.jar
	#sample = /mnt/Storage/home/qinq/softwares/picard-tools-1.91/DownsampleSam.jar ## comment to use built-in sampling tool
	
	[contaminate]    ## for library contamination evaluation, use bowtie for fast assessment only
	
	hg19 = /mnt/Storage/data/Bowtie/hg19
	mm9 = /mnt/Storage/data/Bowtie/mm9
	rn4 = /mnt/Storage/home/qinq/projects/chilin/rn4_index/rn4
	
	[fastqc]
	threads = 5
	
	[bowtie] ## if using bowtie
	max_aign = 1  ## keep unique mappable reads or not
	
	[lib] ## external data
	genome_index = /mnt/Storage/data/Bowtie/hg19                                       ## bowtie or bwa index 
	chrom_bed =  /mnt/Storage/home/qinq/lib/chr_limit_hg19.bed                         ## chromosome limitation BED file
	dhs = /mnt/Storage/data/DHS/DHS_hg19.bed                                           ## your union DHS sites path
	velcro = /mnt/Storage/home/qinq/lib/wgEncodeHg19ConsensusSignalArtifactRegions.bed ## black list
	phast = /mnt/Storage/data/sync_cistrome_lib/conservation/hg19/placentalMammals/    ## Phastcon score from UCSC
	tss = /mnt/Storage/home/qinq/lib/refgenes/hg19.refgene.tss                         ## +- 1kb tss extracted from UCSC refgene, in gcap/static
	
	[hotspot]
	## get from http://www.uwencode.org/proj/hotspot/
	chrom_info = chromosome_info   ## from hotspot website
	mappable_region = path  ## this is downloaded from hotspot website, it's a necessary part for evaluating genomic promotor percentage
	keep_dup = T            ## duplicates or not
	fdrs = "0.01"           ## several fdrs
	
	## if your peak caller is macs2	, fill the following parameters
	[macs2]
	species = hs
	keep_dup = all ## keep all duplicate tags
	shiftsize = 50 ## could be customized

Instructions on the `conf` details.
`Input Format`
support fastq,bam,sam and bed files now.

Fastq Files:

	in the conf files
	`sequence_type` to control sequence type:`se` or `pe`.
	
	If input is single end data, use `,` to separate replicates files.
	If input is pair end data, use `,` to separate pairs, `;` to separate replicates.
	
BAM Files: 
	
	in the conf files, set to `bam, pe` or `bam, se`.
	Query name should be in the neighboring places.
	Original mapping results `SAM` converted by `samtools view -bt` or `picard SortSam SO=queryname` or `samtools sort -n` by query name to make sure that paired reads are in neighboring places, we could use built-in sampling or `[picard] sample` part, because GCAP built-in sampling method only support query name ordered SAM files.
	
	If you are not clear about your mapping parameter, you could try bamToFastq to convert bam to fastq and remapping through our above Fastq scheme.
	
SAM Files, original mapping results with headers , if you only have `bam` files, use `samtools view -h`:

	In the conf file `sequence_type` to `sam, pe` or `sam, se`, files separated by comma.
	
	If you want all SAM files have uniform mapping parameters, you could convert SAM to fastq by samtools view -bt and bamToFastq, then follow up our fastq schemes.

reads BED(converted by bedtools from BAM, sometimes GEO only preserve data with this format):
     
    As our proposals is based on sampling raw reads, including mappable and unmappable reads, BED reads files(BED with 6 fields) do not have unmappable information, so this format is added only for analysis of the rest criteria. As bedToBam could only process SE bed data, PE would be regarded as SE, too. We take all BED format data as SE data, we sample down BED mappable reads 5M for comparison. Change sequence_type to `bed`. This is used for internal data comparison now.
    
    e.g.
    chr1    192388233       192388269       SOLEXA-1GA-2_0072_FC629AV:6:1:3436:1104#0/1     255     -
	chr10   43655355        43655391        SOLEXA-1GA-2_0072_FC629AV:6:1:3567:1104#0/1     255     +
	

`Keep duplicate`
this is an important option for peaks calling. You could customize it by python conf files. To keep duplicates tags,
just `keep_dup = T` in `[hotspot]`.

##### 2. dry run to make sure the installation and conf for needed files are right

* dry run
only print command line: 

	GCAP run -c conf --dry-run
	
If any files are dangling, just refer to install part or email the author.

##### 3. real run 

run with real data: 

	GCAP run -c conf

if any accident happen, debug and resume the program
resume process when problems occurs:

	GCAP run -c conf --resume

* Step control

	- skip steps:

    		GCAP run -c GCAP_pe.conf --skip 9 --resume

	- from and to:
	
		   	GCAP run -c GCAP_pe.conf --from 1 --to 3 --resume
		   	

##### after running successfully
		   	
* clean up and purge
  to remove temporary files, only keep QC report, bigwiggle, peaks, hotspot and bam files, json and latex files:
  
  		GCAP clean -c GCAP_pe.conf
  		GCAP purge -c GCAP_pe.conf
  		
##### batch mode
write a batch file `batch.conf`, which writes the path of the conf files: 
    
    	1.conf
    	2.conf
    	3.conf

the 1.conf, 2.conf, 3.conf is written up to the requirements of above conf, then run:
		
		GCAP batch -b batch.conf --resume

  
#### Tips 
extract Top reads from bam, sam or fastq or Three modes of sampling:

- We use built-in function to do raw reads sampling from PE and SE FASTQ(default) .
- Python function to sample reads from PE and SE SAM(BAM converted SAM) files， including `mappable and unmappable reads`.(default for SAM and BAM)
- picard sampling for PE and SE SAM and BAM files mappable reads. (May use many threads and memory)(optional in picard options, uncomment for default)
- MACS2 sampling for SE BAM files SE mappable reads.(not added)

We decide to use top 5M and 100k mappable and unmappable reads from PE and SE SAM files and both transfer to hospot as single end data.


`samtools view -X` to see the mapping status from column 2nd, see FLAG explanation`http://picard.sourceforge.net/explain-flags.html`
If your SAM/BAM files are not original mapping results, you may need `Restoring pairing information`, this is needed for random access of raw paired reads.

##### sort by name
 samtools sort -n <in.bam> <byname.bam>
##### fix the mate info
 samtools sort fixmate <byname.bam> <byname.fixed.bam>
##### sort by genomic coordinate
 samtools sort <byname.fixed.bam> <out.bam>


### Prototype  Features

1. estimate per sequence quality and library contamination by using 100k sampled reads(mapping by bowtie and bwa(not added yet))
2. reads mapping(exclude mitochrondrial mapping reads) and peaks calling on replicates all reads and 5M sampled reads(peak calling by hotspot and MACS2(added), sampling by Picard DownSampling by probability(5M/total_reads), which seems to be strange, sampling has been replaced by built-in python function, this is an option in conf file. set `picard sample path` would choose picard sampling, other situation would use built-in. Considering sampling from pair 1 for PE fastq or SE fastq files.
3. For pair end data, 5' tags from each pair will be treated as single end for hotspot v3.
4. peaks calling on combo of all reads and 5M sampled readss, Hotspot for 5M reads, Peaks for all reads, that is, use `b, d`.
5. estimate library complexity/redundancy by 5M reads(census;picard, Markduplicates; macs2 filterdup; awk)
6. Add RSC / NSC to QC with SPOT score for 5M reads (Hotspot, a), optional: MACS2 spot score
7. estimate replicates consistency by
	1.a `wigCorrelate` for genome-wide correlation
	1.b `bigWigCorrelate` for top 5M reads BigWiggle Correlation on union DHS regions(filtered by blacklist)
	1.c BigWiggle Correlation on 5M sampled(top 5M) data union hotspot(Hotspot, b, which is filtered to remove blacklist and outlier regions) by bigwiggle(bigWigCorrelate, merged by bedops -m)
    2. Overlap by hotspot(Hotspot, b) regions overlap from 5M reads(Intersection over Union regions, bedtools)
    
8. This has been removed, calculate hotspot(filtered Hotspot, b) promotor percentage and compare with genome promotor percentage for 5M reads Hotspot(b) regions(built-in script)
9. calculate Phastcon score of top 1000 non-promotor Hotspot(filtered Hotspot, b)regions in 100 bp width around summits for 5M reads(modified cistrome built-in module)
10. estimate (narrow peaks, d) overlap with ENCODE narrow peaks union DHS on replicates of 5M reads(bedtools)

Reference
============
1. http://sourceforge.net/apps/mediawiki/srma/index.php?title=User_Guide
2. http://seqanswers.com/forums/showthread.php?t=16375
3. http://picard.sourceforge.net/explain-flags.html


----
Example
=========
- An output from GCAP


![QC report](example.png)
