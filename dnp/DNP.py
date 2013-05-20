#!/usr/bin/env python3

import os
import sys
import re
from glob import glob

import argparse
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from pkg_resources import resource_filename

from dnp.config import Conf
from dnp.helpers import *

latex_template = resource_filename("dnp", "static/Latex_summary_report.jinja2")
## hotspot
token_file = resource_filename("dnp", "static/runall.tokens.txt")
pipeline_scripts = resource_filename("dnp", "pipeline-scripts")


class FriendlyArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit()

def parse_args(args=None):
    """
    If args is None, argparse will parse from sys.argv
    """
    description = "DNP :  DNase I Seq pipeline"
    parser = FriendlyArgumentParser(description=description)
    sub_parsers = parser.add_subparsers(help="sub-command help", dest="sub_command")

    parser_run = sub_parsers.add_parser("run", help="run pipeline using a config file",
        description="DNP-run: Run DNP pipeline using a config file")

    parser_batch = sub_parsers.add_parser("batch", help="run pipeline for multiple datasets")
    parser_batch.add_argument("-b", "--batch-config", dest="batch_config", required=True, help="batch file")
    parser_purge = sub_parsers.add_parser("purge")
    parser_clean = sub_parsers.add_parser("clean", help="Move result file into a new folder and delete other files",
        description="DNP-run: Run DNP pipeline using a config file")

    for p in (parser_run, parser_batch):
        p.add_argument("--from", dest="start_step", default=0, type=int,
            help="Only step after this number will be processed")
        p.add_argument("--to", dest="end_step", default=100, type=int,
            help="Only step before this number will be processed ")
        p.add_argument("--skip", dest="skip_step", default="",
            help="Steps to skip, use comma as seperator")

    for p in (parser_run, parser_batch, parser_purge, parser_clean):
        p.add_argument("-v", "--verbose-level", dest="verbose_level", type=int, default=2)
        p.add_argument("--dry-run", dest="dry_run", action="store_true", default=False)
        p.add_argument("--allow-dangling", dest="allow_dangling", action="store_true", default=False)
        p.add_argument("--resume", dest="resume", action="store_true", default=False)
        p.add_argument("--remove", dest="clean", action="store_true", default=False)

    for p in (parser_run, parser_purge, parser_clean):
        p.add_argument("-c", "--config", dest="config", required=True,
            help="specify the config file to use", )

    return parser.parse_args(args)

def _raw_reads_statistics(workflow, conf):
    """ Kent's codes for raw reads number and reads length
    input: fastq files, PE and SE
    """
    pass

def _sample_reads(workflow, conf, N):
    """
    sample 100k for library contamination
    5M for other evaluation, input fastq
    """
    if N == 100000:
        suffix = "100k"
    elif N == 5000000:
        suffix = "5M"
    if conf.seq_type == "se":
        for n, target in enumerate(conf.treatment_targets):
            attach_back(workflow,
                PythonCommand(single_end_fastq_sampling,
                    input={"fastq": conf.treatment_raws[n]},
                    output={"fastq_sample": target + "_%s.fq" % suffix},
                    param={"random_number": N}))     ## use 100kb random reads
    elif conf.seq_type == "pe":
        for raw, target in conf.treatment_pairs_pe:
            attach_back(workflow,
                PythonCommand(pair_end_fastq_sampling,
                input={"fastq": raw},
                output={"fastq_sample": [ i + "_%s.fq" % suffix for i in target ]},
                param={"random_number": N}))     ## use 100kb random reads

def _fastqc(workflow, conf):
    """
    sample 100k reads for fastqc evaluation
    """
    _sample_reads(workflow, conf, 100000)
    if conf.seq_type == "se":
        for target in conf.treatment_targets:
                fastqc_run = attach_back(workflow,
                    ShellCommand(
                        "{tool} {input[fastq_sample]} --extract -t {param[threads]} -o {output[target_dir]}",
                        input= {"fastq_sample": target + "_100k.fq"},
                        output={"target_dir": conf.target_dir,
                                "fastqc_summary": target + "_100k.fq_fastqc/fastqc_data.txt"},
                        tool="fastqc",
                        param={"threads": 4}))
                fastqc_run.update(param=conf.items("fastqc"))
        ## get per sequence quality median
        attach_back(workflow, PythonCommand(stat_fastqc,
                input = {"fastqc_summaries": [ t + "_100k.fq_fastqc/fastqc_data.txt" for t in conf.treatment_targets ]},
                output = {"json": conf.json_prefix + "_fastqc.json"},
                param = {"ids": conf.treatment_targets}))

    elif conf.seq_type == "pe":
        for n, target in enumerate(conf.treatment_pair_data):
            for p in target:
                fastqc_run = attach_back(workflow,
                    ShellCommand(
                        "{tool} {input[fastq_sample]} --extract -t {param[threads]} -o {output[target_dir]}",
                        input= {"fastq_sample": p + "_100k.fq"},
                        output={"target_dir": conf.target_dir,
                                "fastqc_summary": p + "_100k.fq_fastqc/fastqc_data.txt"},
                        tool="fastqc",
                        param={"threads": 4}))
                fastqc_run.update(param=conf.items("fastqc"))
            attach_back(workflow, PythonCommand(stat_fastqc,
                input = {"fastqc_summaries": [ p + "_100k.fq_fastqc/fastqc_data.txt" for p in target ]},
                output = {"json": conf.json_prefix + "_rep" + str(n+1) + "_fastqc.json"},
                param = {"ids": target}))

def _lib_contamination(workflow, conf):
    """ input PE and SE sampled data 100 K to run,
     output json for library contamination, if not correct, shutdown
    """
    ## bowtie mapping back to mouse, rat and human
    if conf.seq_type == "se":
        for target in conf.treatment_targets:
            for species in dict(conf.items("contaminate")):
                attach_back(workflow,
                    ShellCommand(
                        "{tool} -p {param[threads]} -S -m {param[max_align]} \
                        {param[genome_index]} {input[fastq_sample]} {output[sam]} 2> {output[bowtie_summary]}",
                        input={"fastq_sample": target + "_100k.fq"},
                        output={"sam": target + ".sam.100k",
                                "bowtie_summary": target + species + "_contam_summary.txt"},
                        tool="bowtie",
                        param={"threads": 4,
                               "max_align": 1,
                               "genome_index": conf.get_path("contaminate", species)}))

    elif conf.seq_type == "pe":
        for n, target in enumerate(conf.treatment_pair_data):
            for species in dict(conf.items("contaminate")):
                attach_back(workflow,
                    ShellCommand(
                        "{tool} -X 600 --chunkmbs 300 -n 1 -p {param[threads]} -S -m {param[max_align]} \
                        {param[genome_index]} -1 {input[fastq_sample][0]} -2 {input[fastq_sample][1]} {output[sam]} 2> {output[bowtie_summary]}",
                        input={"fastq_sample": [ t + "_100k.fq" for t in target ]},
                        output={"sam": conf.treatment_targets[n] + ".sam.100k",
                                "bowtie_summary": conf.treatment_targets[n] + species + "_contam_summary.txt"},
                        tool="bowtie",
                        param={"threads": 4,
                               "max_align": 1,
                               "genome_index": conf.get_path("contaminate", species)}))
    all_species = [i for i, _ in conf.items("contaminate")]
    bowtie_summaries = []
    for target in conf.treatment_targets:
        bowtie_summaries.append([target + species + "_contam_summary.txt" for species in all_species])
    attach_back(workflow,
        PythonCommand(stat_contamination,
            input={"bowtie_summaries": bowtie_summaries},
            output={"json": conf.json_prefix + "_contam.json"},
            param={"samples": conf.treatment_bases,
                   "id": conf.id,
                   "species": all_species}))

def _bowtie(workflow, conf):
    """
     all reads mapping
    """
    if conf.seq_type == "se":
        for raw, target in conf.treatment_pairs:
            bowtie = attach_back(workflow,
                            ShellCommand(
                                "{tool} -n 1 -p {param[threads]} -S -m {param[max_align]} \
                                {param[genome_index]} {input[fastq]} {output[sam]} 2> {output[bowtie_summary]}",
                                input={"fastq": raw},
                                output={"sam": target + ".sam.all",
                                        "bowtie_summary": target + "all_um_summary.txt"},
                                tool="bowtie",
                                param={"threads": 4,
                                       "max_align": 1,
                                       "genome_index": conf.get_path("lib", "genome_index")}))
            bowtie.update(param = conf.items("bowtie"))
    elif conf.seq_type == "pe":
        for n, target in enumerate(conf.treatment_targets):
            bowtie = attach_back(workflow,
                        ShellCommand(
                            "{tool} -X 600 --chunkmbs 300 -n 1 -p {param[threads]} -S -m {param[max_align]} \
                            {param[genome_index]} -1 {input[fastq][0]} -2 {input[fastq][1]} {output[sam]} 2> {output[bowtie_summary]}",
                            input={"fastq": conf.treatment_raws[n]},
                            output={"sam": target + ".sam.all",
                                    "bowtie_summary": target + "all_um_summary.txt"},
                            tool="bowtie",
                            param={"threads": 4,
                                   "max_align": 1,
                                   "genome_index": conf.get_path("lib", "genome_index")}))
            bowtie.update(param = conf.items("bowtie"))

    bowtie_summaries = []
    for target in conf.treatment_targets:
        bowtie_summaries.append(target + "all_um_summary.txt")

    attach_back(workflow,
        PythonCommand(stat_bowtie,
            input={"bowtie_summaries": bowtie_summaries},
            output={"json": conf.json_prefix + "_mapping.json"},
            param={"sams": conf.treatment_bases}))



## optionally
def _bwa(workflow, conf):
    """bwa aln /homea6/qinqian/lib/hg19.fa HL1_LNcap_DHT_50U_50_300bp_CTTGTA_L004_R1_001.fastq > HL1_LNcap_DHT_50U_50_300bp_CTTGTA_L004_R1_001.sai"""
    for target in conf.sample_targets:
        attach_back(workflow, ShellCommand(
            "{tool} aln {input[index]} {input[fastq]} > {output[sai]}",
            tool = "bwa",
            input = {"index": conf.get("lib", "genome_index"), "fastq": conf.sample_targets}))

def _macs2(workflow, conf):
    pass

def _reads_mapping(workflow, conf):
    """
    reads mapping for all reads
    and convert sam to bam, sort
    """
    if conf.maptool == "bowtie":
        _bowtie(workflow, conf)
    elif conf.maptool == "bwa":
        _bwa(workflow, conf)
    else:
        print("Not support mapping tool")

    ## SortSam to convert and sort for all reads
    for target in conf.treatment_targets:
        lib_sort = attach_back(workflow, ShellCommand(
            "{tool} -Xmx2g -XX:ParallelGCThreads={param[threads]} -jar {param[sort]} I={input[bam]} O={output[bam]} SO=coordinate \
            VALIDATION_STRINGENCY=SILENT",
            tool = "java",
            input = {"bam": target + ".sam.all"},
            output = {"bam": target + "_picard_sort.bam.all"},
            param = {"threads": 4}))
        lib_sort.update(param = conf.items("picard"))

def _library_complexity(workflow, conf):
    """ use 5M reads to estimate complexity => duplicates
    Gifford's script of Census
    or picard markduplicates
    """
    _sample_reads(workflow, conf, 5000000)

    ## sampling 5M for mapping
    if conf.seq_type == "se":
        for n, target in enumerate(conf.treatment_targets):
            bowtie = attach_back(workflow,
                ShellCommand(
                    "{tool} -n 1 -p {param[threads]} -S -m {param[max_align]} \
                    {param[genome_index]} {input[fastq]} {output[sam]} 2> {output[bowtie_summary]}",
                    input={"fastq": target + "_5M.fq"},
                    output={"sam": target + ".sam.5M",
                            "bowtie_summary": target + "5M_um_summary.txt"},
                    tool="bowtie",
                    param={"threads": 4,
                           "max_align": 1,
                           "genome_index": conf.get_path("lib", "genome_index")}))
            bowtie.update(param = conf.items("bowtie"))
    elif conf.seq_type == "pe":
        for n, target in enumerate(conf.treatment_targets):
            bowtie = attach_back(workflow,
                ShellCommand(
                    "{tool} -X 600 --chunkmbs 300 -n 1 -p {param[threads]} -S -m {param[max_align]} \
                    {param[genome_index]} -1 {input[fastq][0]} -2 {input[fastq][1]} {output[sam]} 2> {output[bowtie_summary]}",
                    input={"fastq":  [ i + "_5M.fq" for i in conf.treatment_pairs_pe[n][1]] },
                    output={"sam": target + ".sam.5M",
                            "bowtie_summary": target + "5M_um_summary.txt"},
                    tool="bowtie",
                    param={"threads": 4,
                           "max_align": 1,
                           "genome_index": conf.get_path("lib", "genome_index")}))
            bowtie.update(param = conf.items("bowtie"))

    # picard markduplicates need SortSam
    for target in conf.treatment_targets:
        lib_sort = attach_back(workflow, ShellCommand(
            "{tool} -Xmx2g -XX:ParallelGCThreads={param[threads]} -jar {param[sort]} I={input[bam]} O={output[bam]} SO=coordinate \
            VALIDATION_STRINGENCY=SILENT",
            tool = "java",
            input = {"bam": target + ".sam.5M"},
            output = {"bam": target + "_picard_sort.bam.5M"},
            param = {"threads": 4}))

        lib_sort.update(param = conf.items("picard"))
    # markduplicates
    for target in conf.treatment_targets:
        lib_dup = attach_back(workflow, ShellCommand(
            "{tool} -Xmx2g -XX:ParallelGCThreads={param[threads]} -jar {param[markdup]} I={input[bam]} O={output[bam]} METRICS_FILE={output[metrics]} REMOVE_DUPLICATES=false \
            VALIDATION_STRINGENCY=SILENT",
            tool = "java",
            input = {"bam": target + "_picard_sort.bam.5M"},
            output = {"bam": target + "_markdup.bam.5M", "metrics": target + "_markdup_metric.5M"},
            param = {"threads": 4}))
        lib_dup.update(param = conf.items("picard"))

def _fragment_size(workflow, conf):
    """
    Pair end picard, Single end macs2
    """
    if conf.seq_type == "pe":
        for target in conf.treatment_targets:
            insert_size = attach_back(workflow, ShellCommand(
                        "{tool} -Xmx2g -XX:ParallelGCThreads={param[threads]} -jar {param[insertsize]} HISTOGRAM_FILE={output[pdf]} I={input[bam]} O={output[insert]} VALIDATION_STRINGENCY=SILENT",
                        tool = "java",
                        input = {"bam": target + "_picard_sort.bam.5M"},
                        output = {"pdf": target + "_picard_insert_5M.pdf", "insert": target + "_insert_metric.5M"},
                        param = {"threads": 4}))
            insert_size.update(param = conf.items("picard"))
    elif conf.seq_type == "se":
        for target in conf.treatment_targets:
            fragment_size = attach_back(workflow, ShellCommand(
                        "{tool} predictd -i {input[bam]} --rfile {param[prefix]}",
                        tool = "macs2",
                        input = {"bam": target + "_picard_sort.bam.5M"},
                        output = {"R": target + "_5M_model.R"},
                        param = {"threads": 4, "prefix": target + "_5M"}))
            fragment_size.update(param = conf.items("picard"))

def _hotspot_on_replicates(workflow, conf):
    # I write a tags.sh in pipeline-scripts to convert BAM files to 1 bp BED files
    # for both SE and PE
    # using all and 5M reads
    ## all reads, need to get 1bp tags density, and 20bp tags density to bigwiggle
    ## use two passes merge peaks and hotspot to estimate peaks number
    ## sometimes wavepeaks don't get peaks evenly distributed, shut down chromosome check
    for i, target in enumerate(conf.treatment_targets):
        ## prepare input for hotspot, get mappable tags and starch into output directory
        attach_back(workflow,
            ShellCommand(
                "{tool} {input} {output} {param[type]}",
                tool = "tags.sh",
                input = target + "_picard_sort.bam.all",
                output = target + ".bed.starch",
                param = {"type": conf.seq_type}))
        ## generate configuration for hotspot v3
        hotspot_conf = attach_back(workflow,
            PythonCommand(spot_conf,
                input = {"spot_conf": token_file,
                         "tag": target + ".bed.starch",
                         "mappable_region": conf.get("hotspot", "mappable_region"),
                         "chrom_info": conf.get("hotspot", "chrom_info")},
                output = {"conf": target + "_runall.tokens.txt",
                          "dir": conf.target_dir},
                param = {"fdrs": "0.01", "K": conf.get("Basis", "read_length"),
                         "species": conf.get("Basis", "species")}))
        hotspot_conf.param.update(conf.items("hotspot"))
        ## run hotspot v3
        attach_back(workflow, ShellCommand(
            "{tool} {input[pipe]} {input[token]} {output[dir]}",
            tool = "runhotspot",
            input = {"pipe": pipeline_scripts, "token": target + "_runall.tokens.txt"},
            output = {"dir": conf.target_dir,
                      "spot_peaks_combined": conf.hotspot_reps_both_passes_prefix[i] + ".hotspot.twopass.fdr0.01.merge.pks.bed",
                      "density_starch": target + ".tagdensity.bed.starch",
                      "spot": target + ".spot.out"},
            param = None))

        ## TODO: ## pileup with macs2 extsize 1, using 1bp target + ".bed.starch.all"
        ## for cutting bias and genome browser

    ## all reads, calculate the merged two passes and peaks regions number
    ## Use 5M reads to estimate SPOT
    for i, target in enumerate(conf.treatment_targets):
        ## prepare input for hotspot, get mappable tags and starch into output directory
        attach_back(workflow,
            ShellCommand(
                "{tool} {input} {output} {param[type]}",
                tool = "tags.sh",
                input = target + "_picard_sort.bam.5M",
                output = target + "_5M.bed.starch",
                param = {"type": conf.seq_type}))
        ## generate configuration for hotspot v3
        hotspot_conf = attach_back(workflow,
            PythonCommand(spot_conf,
                input = {"spot_conf": token_file,
                         "tag": target + "_5M.bed.starch",
                         "mappable_region": conf.get("hotspot", "mappable_region"),
                         "chrom_info": conf.get("hotspot", "chrom_info")},
                output = {"conf": target + "_runall_5M.tokens.txt",
                          "dir": conf.target_dir},
                param = {"fdrs": "0.01", "K": conf.get("Basis", "read_length"),
                         "species": conf.get("Basis", "species")}))
        hotspot_conf.param.update(conf.items("hotspot"))
        ## run hotspot v3
        attach_back(workflow, ShellCommand(
            "{tool} {input[pipe]} {input[token]} {output[dir]}",
            tool = "runhotspot",
            input = {"pipe": pipeline_scripts, "token": target + "_runall_5M.tokens.txt"},
            output = {"dir": conf.target_dir,
                      "spot_peaks_combined": conf.hotspot_reps_both_passes_5M_prefix[i] + ".hotspot.twopass.fdr0.01.merge.pks.bed",
                      "density_starch": target + "_5M.tagdensity.bed.starch",
                      "spot": target + "_5M.spot.out"},
            param = None))

        ## density from hotspot 3, 20bp resolution
        ## for correlation evaluation
        attach_back(workflow, ShellCommand(
            "{tool} {input[starch]} > {output[bed]}",
            tool = "unstarch",
            input = {"starch": target + "_5M.tagdensity.bed.starch"},
            output = {"bed": target + "_5M.tagdensity.bed.tmp"}))

        ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
        ## filter bdg file to remove over-border coordinates
        ## with bedtools or bedClip
        attach_back(workflow,
            ShellCommand(
                '{tool} intersect -a {input} -b {param[chrom_bed]} -f 1.00 > {output}',
                tool="bedtools",
                input=target + "_5M.tagdensity.bed.tmp",
                output=target + "_5M.tagdensity.bed",
                param={'chrom_bed': conf.get("lib", "chrom_bed")},
                name="bed replicate filtering"))
        ## convert to bigwiggle
        attach_back(workflow,
            ShellCommand(
                'cut -f 1,2,3,5 {input} > {input}.tmp && {tool} {input}.tmp {param[chrom_len]} {output}',
                tool = "bedGraphToBigWig",
                input=target + "_5M.tagdensity.bed",
                output=target + "_5M.bw",
                param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))
    ## 5M reads, calculate the merged two passes and peaks regions number
    attach_back(workflow, PythonCommand(
        stat_spot_on_replicates,
        input = {"spot_files": [ target + "_5M.spot.out" for target in conf.treatment_targets ]},
        output = {"json": conf.json_prefix + "_sample_spot_5M.json"}, ## try density plot
        param = None))

def _hotspot_reps_preprocessing(workflow, conf):
    """
    remove blacklist and outlier
    """
    for i, target in enumerate(conf.treatment_targets):
        ## remove blacklist
        attach_back(workflow,
            ShellCommand(
                "{tool} -v -a {input[peaks]} -b {input[velcro_peaks_bed]} > {output}",
                tool="intersectBed",
                input={"peaks": conf.hotspot_reps_both_passes_5M_prefix[i] + ".hotspot.twopass.fdr0.01.merge.pks.bed",
                       "velcro_peaks_bed": conf.get("lib", "velcro")},
                output = target + "_5M_velcro_non_overlap_peaks.bed",
                name = "filter blacklist",
                param=None))
        ## awk and bedClip to remove outlier location from above input
        attach_back(workflow,
            ShellCommand(
                "{tool} '{{if ($2 >= 0 && $2 < $3) print}}' {input} > {output}",
                tool="awk",
                input = target + "_5M_velcro_non_overlap_peaks.bed",
                output = target + "_5M_velcro_non_overlap_peaks.bed.tmp",
                name = "filter bed files outlier location"))
        # prototype used here to do the similar thing on bedclip
        bed_clip = attach_back(workflow,
            ShellCommand(
                template="{tool} {input} {param[chrom_len]} {output}",
                tool="bedClip",
                input=target + "_5M_velcro_non_overlap_peaks.bed.tmp",
                output = target + "_5M_velcro_non_overlap_peaks.bed.final",
                param={'chrom_len': conf.get_path("lib", "chrom_len")},
                name="bedclip filter"))
        bed_clip.allow_fail = True

def _hotspot_reps_evaluation(workflow, conf):
    ## use 5M reads to evaluate
    ## preprocessing

    ## intersect region ratio in pairwise ways
    for i in range(len(conf.treatment_targets)):
        for j in range(i+1, len(conf.treatment_targets)):
            attach_back(workflow,
                ShellCommand(
                    "{tool} -a {input[one_rep]} -b {input[another_rep]} > {output}",
                    tool="intersectBed",
                    input = {"one_rep": conf.treatment_targets[i] + "_5M_velcro_non_overlap_peaks.bed.final",
                             "another_rep": conf.treatment_targets[j] + "_5M_velcro_non_overlap_peaks.bed.final"},
                    output= conf.treatment_targets[i] + "_reps_" + str(i) + str(j) + "_5M_velcro_non_overlap_peaks.bed.final",
                    name = "replicates Overlap",
                    param=None))

    ## correlation in union regions
    ## union peaks regions
    union_peaks = attach_back(workflow,
        ShellCommand(
            "{tool} -m {param[rep_peaks]} > {output}",
            tool = "bedops",
            input = [ target + "_5M_velcro_non_overlap_peaks.bed.final" for target in conf.treatment_targets ],
            output = conf.prefix + "_reps_union_region.bed",
            name = "union replicates peaks",
            param = {"rep_peaks": " ".join([ target + "_5M_velcro_non_overlap_peaks.bed.final" for target in conf.treatment_targets ])}))
    cor_on_bw = attach_back(workflow,
        ShellCommand(
            """{tool} --min-score {param[wig_correlation_min]} --max-score {param[wig_correlation_max]} \
            -r {output[R]} {param[bw]} -b {input[bed]} {param[rep]} 1> {output[cor]}""",
            tool="bigwig_correlation_in_bed_file.py",
            input= {"bw": [ target + "_5M.bw" for target in conf.treatment_targets ],
                    "bed": conf.prefix + "_reps_union_region.bed"},
            output={"R": conf.prefix + "_cor.R",
                    "cor": conf.prefix + "_5M_cor_stdout"},
            param={"wig_correlation_method": "mean",
                   "wig_correlation_min": 2,
                   "wig_correlation_max": 50,
                   "wig_correlation_step": 10},
            name = "Correlation on union replicates peaks"))
    cor_on_bw.param["bw"] = " ".join([ " -w " +  i for i in cor_on_bw.input["bw"]])
    cor_on_bw.param["rep"] = " ".join([" -l replicate_%s" % (x + 1) for x in range(len(conf.treatment_pairs))])
    cor_on_bw.update(param=conf.items("correlation"))
    cor_on_bw.allow_fail = True
    attach_back(workflow,
        PythonCommand(
            stat_cor,
            input={"correlation_R": conf.prefix + "_cor.R",
                   "cor": conf.prefix + "_5M_cor_stdout"},
            output={"json": conf.json_prefix + "_cor.json"}))

def _hotspot_combo(workflow, conf):
    ## starchcat all replicates
    ## all reads combo
    attach_back(workflow,
        ShellCommand(
            "{tool} {param[tag_list]} {output}",
            tool = "starchcat",
            input = [ target + ".bed.starch" for target in conf.treatment_targets ],
            output = conf.prefix + "_merge_all.bed.starch",
            param = {"tag_list": " ".join([ target + ".bed.starch" for target in conf.treatment_targets ])}))

    ## config for hotspot 3
    hotspot_conf = attach_back(workflow,
        PythonCommand(spot_conf,
            input = {"spot_conf": token_file,
                     "tag": conf.prefix + "_merge_all.bed.starch",
                     "mappable_region": conf.get("hotspot", "mappable_region"),
                     "chrom_info": conf.get("hotspot", "chrom_info")},
            output = {"conf": conf.prefix + "_runall_all.tokens.txt",
                      "dir": conf.target_dir},
            param = {"fdrs": "0.01", "K": conf.get("Basis", "read_length"),
                     "species": conf.get("Basis", "species")}))
    hotspot_conf.param.update(conf.items("hotspot"))

    ## run hotspot v3
    attach_back(workflow, ShellCommand(
        "{tool} {input[pipe]} {input[token]} {output[dir]}",
        tool = "runhotspot",
        input = {"pipe": pipeline_scripts, "token": conf.prefix + "_runall_all.tokens.txt"},
        output = {"dir": conf.target_dir,
                  "spot_peaks_combined": conf.hotspot_merge_both_passes_prefix + ".hotspot.twopass.fdr0.01.merge.pks.bed",
                  "density_starch": conf.prefix + ".tagdensity.bed.starch",
                  "spot": conf.prefix + ".spot.out"},
        param = None))

    ## 5M read combo
    starch_cat = attach_back(workflow,
                    ShellCommand(
                        "{tool} {param[tag_list]} {output}",
                        tool = "starchcat",
                        input = [ target + "_5M.bed.starch" for target in conf.treatment_targets ],
                        output = conf.prefix + "_merge_5M.bed.starch",
                        param = {}))
    starch_cat.param["tag_list"] = " ".join(starch_cat.input)

    ## config for hotspot 3
    hotspot_conf = attach_back(workflow,
        PythonCommand(spot_conf,
            input = {"spot_conf": token_file,
                     "tag": conf.prefix + "_merge_5M.bed.starch",
                     "mappable_region": conf.get("hotspot", "mappable_region"),
                     "chrom_info": conf.get("hotspot", "chrom_info")},
            output = {"conf": conf.prefix + "_runall_5M.tokens.txt",
                      "dir": conf.target_dir},
            param = {"fdrs": "0.01", "K": conf.get("Basis", "read_length"),
                     "species": conf.get("Basis", "species")}))
    hotspot_conf.param.update(conf.items("hotspot"))

    ## run hotspot v3
    attach_back(workflow, ShellCommand(
        "{tool} {input[pipe]} {input[token]} {output[dir]}",
        tool = "runhotspot",
        input = {"pipe": pipeline_scripts, "token": conf.prefix + "_runall_5M.tokens.txt"},
        output = {"dir": conf.target_dir,
                  "spot_peaks_combined": conf.hotspot_merge_5M_both_passes_prefix + ".hotspot.twopass.fdr0.01.merge.pks.bed",
                  "density_starch": conf.prefix + "_5M.tagdensity.bed.starch",
                  "spot": conf.prefix + "_5M.spot.out"},
        param = None))

def _promoter_overlap(workflow, conf):
    """ use bedtools to get promotor peaks enrichment
    use 5M reads peaks region to get this value
    use CEAS with BED file to estimate ChIP promotor proportion
    """
    for i, target in enumerate(conf.treatment_targets):
        attach_back(workflow,
            ShellCommand(
                "{tool} -b {input[bed]} -g {input[gene_table]} --name {param[name]}",
                tool = "ceas",
                input = {"bed": target + "_5M_velcro_non_overlap_peaks.bed.final",
                         "gene_table": conf.get("lib", "gene_table")},
                output = {"R": target + "_5M_velcro_non_overlap_peaks.R",
                          "pdf": target + "_5M_velcro_non_overlap_peaks.pdf"},
                param = {"name": target + "_5M_velcro_non_overlap_peaks"},
                name = "promotor overlaps"))
        ## extract CEAS
        attach_back(workflow, PythonCommand(
            stat_promotor,
            input = {"ceas_r": target + "_5M_velcro_non_overlap_peaks.R"},
            output = {"json": conf.json_prefix + "_rep_" + str(i) + "_promotor.json"}))
    if len(conf.treatment_pairs) >= 2:
        attach_back(workflow,
            ShellCommand(
                "{tool} -b {input[bed]} -g {input[gene_table]} --name {param[name]}",
                tool = "ceas",
                input = {"bed": conf.hotspot_merge_5M_both_passes_prefix + ".hotspot.twopass.fdr0.01.merge.pks.bed",
                         "gene_table": conf.get("lib", "gene_table")},
                output = {"R": target + "_5M_velcro_non_overlap_peaks.R",
                          "pdf": target + "_5M_velcro_non_overlap_peaks.pdf"},
                param = {"name": target + "_5M_velcro_non_overlap_peaks"},
                name = "promotor overlaps"))

def _conservation(workflow, conf):
    ## use 5M reads for evaluation
    for target in conf.treatment_targets:
        ## get non-promotor peaks regions
        attach_back(workflow,
            ShellCommand(
                "{tool} -v -a {input[peaks]} -b {input[tss]} > {output}",
                tool = "intersectBed",
                input = {"peaks": target + "_5M_velcro_non_overlap_peaks.bed.final", "tss": conf.get("lib", "tss")},
                output = target + "_5M_non-promotor_peaks_bed",
                name = "Get non-promotor peaks regions",
                param = None))

        get_top_peaks = attach_back(workflow,
            ShellCommand(
                "{tool} -r -g -k 5 {input} | head -n {param[peaks]} > {output}",
                tool="sort",
                input=target + "_5M_non-promotor_peaks_bed",
                output=target + "_5M_non_promotor.bed.top1000",
                param={'peaks': 1000}, name="top summits for conservation"))
        get_top_peaks.update(param=conf.items('conservation'))

        ## get summits is implemented in conservation_average.py
        attach_back(workflow,
            ShellCommand(
                "{tool} -w {param[width]} -d {input[phastcon_db]} {input[bed]} 1>{output}",
                tool = "conservation_average.py",
                input= {"phastcon_db": conf.get_path("lib", "phast"),
                        "bed": target + "_5M_non_promotor.bed.top1000"},
                output = target + "_100bp_phastcon.score",
                param = {"width": 100}))

        attach_back(workflow,
            PythonCommand(stat_conserv,
                input = target + "_100bp_phastcon.score",
                output = {"json": target + "_top1000_100bp_phast.json"},
                param = None))

def _union_DHS_overlap(workflow, conf):
    for target in conf.treatment_targets:
        attach_back(workflow,
            ShellCommand(
                "{tool} -wa -u  \
                -a {input[pks_spot_bed]} -b {input[DHS_peaks_bed]} > {output}",
                tool="intersectBed",
                input={"pks_spot_bed": target + "_5M_velcro_non_overlap_peaks.bed.final",
                       "DHS_peaks_bed": conf.get("lib", "dhs")},
                output=conf.prefix + "_DHS_overlap_peaks_bed",
                name = "Write out DHS overlap BED",
                param=None))
        attach_back(workflow, PythonCommand(
            stat_dhs,
            input={"dhs_peaks": conf.prefix + "_DHS_overlap_peaks_bed",
                   "pks_spot_bed":  target + "_5M_velcro_non_overlap_peaks.bed.final"},
            output={"json": conf.json_prefix + "_dhs.json"},
            param=None,
            name="DHS summary"))

def prepare_clean_up(workflow, conf):
    """
    package all the necessary results and delete temporary files
    preserve bed, bigwiggle, starch, bam and pdf report
    """
    p_list = ['*.bam', '*.xls', '*_summits.bed', '*_peaks.bed', '*.bw',
              '*.R', '*.zip', '*cor*', 'json', "*summary*",
              "*seqpos", "*fastqc", '*latex', "*.conf"]

    p_pattern = [os.path.join(conf.target_dir, p) for p in p_list]

    final_dir = conf.target_dir + '/dataset_' + conf.id
    attach_back(workflow,
        ShellCommand("if [ ! -d '{output}' ]; then mkdir -p {output}; fi",
            output=final_dir))

    for pf in p_pattern:
        if not glob(pf):
            print(pf)
            continue
        move = attach_back(workflow,
            ShellCommand('mv {param[preserve_files]} {output[dir]} \n# Pattern: {param[p_pattern]}',
                output={"dir": final_dir},
                param={"preserve_files": " ".join(glob(pf)),
                       "p_pattern": pf}, ))
        move.allow_fail = True

def prepare_purge(workflow, conf):
    d_list = ["json", "latex", conf.id + "_whole.*"]
    d_pattern = [os.path.join(conf.target_dir, d) for d in d_list]

    for df in d_pattern:
        if not glob(df):
            print(df)
            continue
        deleted = attach_back(workflow,
            ShellCommand("rm -r {param[deleted_files]}",
                param={"deleted_files": " ".join(glob(df))}))
        deleted.allow_fail = True

def create_workflow(args, conf):
    workflow = Workflow(name="Main")
    bld = DNPBuilder(workflow, conf)

    if args.sub_command == "clean":
        ## not implemented yet
        bld.build(prepare_clean_up)
        return workflow

    if args.sub_command == "purge":
        bld.build(prepare_purge)
        return workflow

    if args.skip_step:
        skipped_steps = [int(i) for i in args.skip_step.split(",")]
    else:
        skipped_steps = []

    step_checker = StepChecker(args.start_step, args.end_step, skipped_steps)

    have_treat_reps = len(conf.treatment_pairs) >= 2 ## replicates
    has_dhs = conf.get("lib", "dhs")
    has_velcro = conf.get("lib", "velcro")  ## mouse not have blacklist yet
    need_run = step_checker.need_run

    bld.attach_back(ShellCommand(
        "if [ ! -d '{output}' ]; then mkdir -p {output}; fi",
        output=conf.target_dir))

    if need_run(1):
        bld.build(_fastqc)

    if need_run(2):
        bld.build(_lib_contamination)

    if need_run(3):
        bld.build(_reads_mapping)

    if need_run(4):
        bld.build(_library_complexity)

    if need_run(5):
        bld.build(_fragment_size)

    if need_run(6):
        bld.build(_hotspot_on_replicates)
        bld.build(_hotspot_reps_preprocessing)
        if have_treat_reps:
            bld.build(_hotspot_reps_evaluation)
            bld.build(_hotspot_combo)
    if need_run(7):
        bld.build(_promoter_overlap)
    if need_run(8):
        bld.build(_conservation)
    if has_dhs and need_run(9):
        bld.build(_union_DHS_overlap)
#    if need_run(10):
#        bld.build(_report())
    return workflow

def main(args=None):
    args = parse_args(args)
    print("Arguments:", args)

    if args.sub_command in ["run", "clean", "purge"]:
        conf = Conf(args.config)
        workflow = create_workflow(args, conf)

        workflow.set_option(
            verbose_level=args.verbose_level,
            dry_run_mode=args.dry_run,
            resume=args.resume,
            allow_dangling=args.allow_dangling)

        workflow.invoke()

    elif args.sub_command == "batch":
        with open(args.batch_config) as batch_file:
            for a_conf in batch_file:
                a_conf = a_conf.strip()
                conf = Conf(a_conf)
                workflow = create_workflow(args, conf)
                workflow.set_option(
                    verbose_level=args.verbose_level,
                    dry_run_mode=args.dry_run,
                    resume=args.resume,
                    allow_dangling=args.allow_dangling)
                workflow.invoke()

if __name__ == "__main__":
    main()
