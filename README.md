## What is DNP
> ENCODE 3 DNase I seq QC pipeline

Install 
=============
hg clone https://bitbucket.org/hanfeisun/samflow
or
pip install samflow

install bowtie
----------------
add binary to PATH and download hg19, mm9, rn4 index for library contamination evaluation

install hotspot
----------------
Download hotspot v3 from , get the hotspot5, wavelets, wavePeaks into `PATH`

install bedops and bedtools
-----------------------------

***

install MACS2
------------------
install as follows: 
	
	git clone https://github.com/taoliu/MACS/
	or 
	pip install MACS2

***

install picard
--------------------
Use picard for SortSam, Markduplicates for both single end and pair end data.
Use picard for pair end data `fragment size` evaluation.
For single end data, we used `MACS2 predictd`

----

install cistrome-application
----------------------------
Use `mercurial` to get necessary packages:

    hg clone https://bitbucket.org/cistrome/cistrome-applications-harvard
    
install each of the packages contained in `cistrome-applications-harvard` repo

* download gene_tables for CEAS and Phastcon score bigwiggles of Placetalmammals for human and mouse
* use the lastest repo

- - - -

Usage
=============

refer to static/DNP_pe.conf for pair end data, static/DNP_se.conf for single end data.

If input is single end data, use `,` to separate replicates files.
If input is pair end data, use `,` to separate pairs, `;` to separate replicates.

instructions on conf files:

	[lib] 
	tss = 
for human should be the data contained in the static/hg19.refgene.tss


* dry run

only print command line: 

	DNP.py run -c DNP_se.conf

* real run

run with real data: 

	DNP.py run -c DNP_pe.conf

* resume

resume process when problems occurs:

	DNP.py run -c DNP_pe.conf --resume


* skip steps
DNP.py run -c DNP_pe.conf --skip 9 --resume

* from and to

DNP.py run -c DNP_pe.conf --from 1 --to 3 --resume

----

Prototype  Features
======================
Following features has been tested in single end data for single replicates successfully.

1. estimate library contamination by using 100k sampled reads

2. peaks calling on replicates all reads and 5M sampled reads

3. peaks calling on combo of all reads and 5M sampled reads

3. estimates SPOT score for 5M reads peaks

4. estimate replicates consistency by using 5M sampled data union peaks by bigwiggle and by overlaps

5. calculate top 1000 non-promotor peaks regions in 100 bp width around summits.

6. estimate peaks overlap with DHS on replicates and combo by 5M reads

TODO
=====================
1. Add bwa for reads mapping optionally
2. Add macs2 for peaks calling optionally
3. add a jinja2 template for QC summary latex report rendering, e.g. for pair end, 
4. test on more data and situations, like pair end data with replicates.

![QC report](example.png)
