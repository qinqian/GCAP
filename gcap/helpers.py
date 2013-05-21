from configparser import ConfigParser
import json
import os
import random
import re
import tempfile
from jinja2 import Environment, FileSystemLoader
import math
from samflow.command import AbstractCommand, ShellCommand
from samflow.workflow import attach_back
import sys

class StepChecker:
    def __init__(self, start, end, skips):
        self.start = start
        self.end = end
        self.skips = skips

    def need_run(self, step_id):
        if step_id < self.start:
            return False
        if step_id > self.end:
            return False
        if step_id in self.skips:
            return False
        return True

class GCAPBuilder:
    def __init__(self, workflow, conf):
        self.workflow = workflow
        self.conf = conf
        self.LaTex_fragments = []
        self.plain_fragments = []
        self.finished = set()

    def build(self, prep_func, tag=None):
        prep_func(self.workflow, self.conf)
        if tag:
            self.finished.add(tag)

    def attach_back(self, command):
        attach_back(self.workflow, command)

def json_dump(json_dict):
    json_file = json_dict["output"]["json"]
    with open(json_file, "w") as f:
        json.dump(json_dict, f, indent=4)
    return json_file

def json_load(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    return json_dict

env = Environment(loader = FileSystemLoader("/"),
    block_start_string = '\BLOCK{',
    block_end_string = '}',
    variable_start_string = '\VAR{',
    variable_end_string = '}',
    comment_start_string = '\#{',
    comment_end_string = '}',
    line_statement_prefix = '%-',
    line_comment_prefix = '%#',
    trim_blocks = True,
    autoescape = False,)

def surround_by_quote(a_list):
    return ['"%s"' % an_element for an_element in a_list]

def count_in_million(x):
    if type(x) == int:
        return str(round(x/1000000., 3)) + "M"

def decimal_to_latex_percent(dec):
    if type(dec) == float:
        return str(round(dec/100, 3))    + "\%"

def underline_to_space(x):
    if type(x) == str:
        return x.replace("_", " ")
    return x

env.filters["surround_by_quote"] = surround_by_quote

class JinjaTemplateCommand(AbstractCommand):
    def __init__(self, template, tool=None, param = {}, input=[], output=[], name = ""):
        AbstractCommand.__init__(self, template=template, tool=None, param = param, input=[], output=[], name = name)
        self.env = env
        self._t = self.env.get_template(self.template)

    def _execute(self):
        """ Real-run current command"""
        self.result = self._t.render(input = self.input, output = self.output, **self.param)
        return True

    def _simulate(self):
        """ Dry-run current command: Pretend to run but not invoke anything """
        print("Rendering Latex part %s" % self.name, self.template)
        return True

def template_dump(jinja_template):
    jinja_template.invoke()
    with open(jinja_template.param["render_dump"], "w") as f:
        f.write(jinja_template.result)

def single_end_fastq_sampling(input = {"fastq": ""}, output = {"fastq_sample": ""}, param = {"random_number": ""}):
    num_lines = sum(1 for _ in open(input["fastq"])) / 4

    rand_nums = sorted([random.randint(0, num_lines - 1) for _ in range(param["random_number"])])

    fastq = open(input["fastq"],"rU")
    fastq_sample = open(output["fastq_sample"], "w")

    cur_num = - 1
    for rand_num in rand_nums:
        while cur_num < rand_num:
            for i in range(4):
                fastq.readline()
            cur_num += 1

        for i in range(4):
            fastq_sample.write(fastq.readline())
        cur_num += 1

    fastq_sample.close()
    fastq.close()


def pair_end_fastq_sampling(input = {"fastq": ""}, output = {"fastq_sample": ""}, param = {"random_number": ""}):
    def write_random_records(fqa, fqb, suba, subb, N=100000):
        records = sum(1 for _ in open(fqa)) / 4
        rand_records = sorted([random.randint(0, records - 1) for _ in range(N)])

        fha, fhb = open(fqa),  open(fqb)
        suba, subb = open(suba, "w"), open(subb, "w")
        rec_no = - 1
        for rr in rand_records:
            while rec_no < rr:
                rec_no += 1
                for i in range(4): fha.readline()
                for i in range(4): fhb.readline()
            for i in range(4):
                suba.write(fha.readline())
                subb.write(fhb.readline())
            rec_no += 1
        suba.close()
        subb.close()
        fha.close()
        fhb.close()
    write_random_records(input["fastq"][0], input["fastq"][1], output["fastq_sample"][0], output["fastq_sample"][1], param["random_number"])

## summary of library contamination
def stat_contamination(input = {"bowtie_summaries": [[]]},
                       output = {"json": ""}, param = {"samples": "", "species": "", "id": ""}):

    library_contamination = {}
    library_contamination["header"] = {"sample": "sample name", "species": param["species"]}
    library_contamination["value"] = {}
    for a_summary, s in zip(input["bowtie_summaries"], map(underline_to_space, param["samples"])):
        ## each bowtie_summary has several species information
        library_contamination["value"][s] = {}
        for i, j in zip(a_summary, param["species"]):
            ## species 1, species2, species3
            parsed_summary = _bowtie_summary_parse(i)
            library_contamination["value"][s][j] = parsed_summary["mappable_rate"]

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict["stat"] = library_contamination
    json_dump(json_dict)

def fastqc_parse(input):
    data = open(input).readlines()
    sequence_length = 0
    quality_dict = {}
    in_seq_quality_section = False
    for line in data:
        if re.search(r"^Sequence length", line):
            assert sequence_length == 0
            sequence_length = int(re.findall(r"^Sequence length\t(\d+)", line)[0])
        elif re.search(r"^>>Per sequence quality", line):
            assert not in_seq_quality_section
            in_seq_quality_section = True
            continue

        if re.search(r"^>>END_MODULE", line) and in_seq_quality_section:
            in_seq_quality_section = False

        if (not line.startswith("#")) and in_seq_quality_section:
            sequence_quality = re.findall("^(\w+)\t(\w+)", line)[0]
            quality_dict[sequence_quality[0]] = float(sequence_quality[1])
    total = sum(quality_dict.values())
    n = 0
    for item in sorted(quality_dict.items(), key=lambda e: e[0], reverse=True):
        n = n + item[1]
        if n / total > 0.5:
            median = int(item[0])
            break
    return {"sequence_length": sequence_length,
            "median": median}

def stat_fastqc(input = {"fastqc_summaries": []},
                output={"json": ""},
                param= {"ids": ""}):
    """
    per sequence quality
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    stat = {}
    for a_summary, a_id in zip(input["fastqc_summaries"], param["ids"]):
        parsed = fastqc_parse(input=a_summary)
        stat[a_id] = {}
        stat[a_id]["median"] = parsed["median"]
        stat[a_id]["cutoff"] = 25
        stat[a_id]["sequence_length"] = parsed["sequence_length"]
    print(stat)
    json_dict["stat"] = stat
    json_dump(json_dict)

def _bowtie_summary_parse(input=""):
    summary_content = open(input).read()
    print(summary_content)
    print("_" * 100)
    total_reads = int(re.findall(r"# reads processed: (.*)\n", summary_content)[0])

    # WARN: this is `unique mappable reads` as we set `-m 1`
    mappable_reads = int(re.findall(r"# reads with at least one reported alignment: (.*) \(", summary_content)[0])

    mappable_rate = float(mappable_reads) / float(total_reads)

    return {"total_reads": total_reads, "mappable_reads": mappable_reads, "mappable_rate": mappable_rate}

def stat_bowtie(input={"bowtie_summaries": []},
                output={"json": ""},
                param={"sams": []}):
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}

    for summary, sam in zip(input["bowtie_summaries"], param["sams"]):
        json_dict["stat"][sam] = _bowtie_summary_parse(summary)
        json_dict["stat"][sam]["cutoff"] = 0.5 # mappable reads
    json_dump(json_dict)


def stat_markdup(input = {"markdup": ""},
                 output = {"json": ""},
                 param = None):
    data = open(input["markdup"]).readlines()
    dup = data[7].split("\t")
    json_dict = {"input": input, "output": output, "param": param, "stat": dup}
    json_dump(json_dict)

def spot_conf(input = {"tag": "", "mappable_region": "", "spot_conf": "", "chrom_info": ""},
              output = {"dir": "", "conf": ""},
              param = {"K": "", "FDRS": "0.01", "species": ""}):
    cf = ConfigParser()
    ## preserve uppercase
    cf.optionxform = str
    cf.read(input["spot_conf"])
    setting = lambda cf, before, after: cf.set("script-tokenizer", before, after)
    setting(cf,"_TAGS_", input["tag"])
    setting(cf,"_GENOME_", param["species"])
    setting(cf, "_K_", param["K"])
    setting(cf,"_MAPPABLE_FILE_", input["mappable_region"])
    setting(cf, "_CHROM_FILE_", input["chrom_info"])
    setting(cf, "_DENS_", output["dir"])
    setting(cf, "_OUTDIR_", output["dir"])
    setting(cf, "_RANDIR_", output["dir"])
    setting(cf, "_FDRS_", param["fdrs"])
    setting(cf, "_DUPOK_", "T")
    setting(cf, "_CHKCHR_", "chrX")     ## To check if data contains this chromosome
    setting(cf, "_CHECK_", "F")         ## change to F in case that hotspot error for no wavelets peaks, other method
    cf.write(open(output["conf"], "w"))

def stat_spot(input = {"spot_files": ""}, output = {"json": ""},
              param = None):
    spot = {}
    spot["input"] = input
    spot["output"] = output
    print(input["spot_files"])
    spot["stat"] = {}
    with open(input["spot_files"], "rU") as f:
        stat = f.readlines()[1].strip().split()
        spot["stat"]["total"] = float(stat[0])
        spot["stat"]["hotspot"] = float(stat[1])
        spot["stat"]["spot"] = float(stat[2])
    json_dump(spot)
    return

def stat_spot_on_replicates(input = {"spot_files": ""}, output = {"json": "", "R": "", "pdf": ""},
                            param = None):
    spot = {}
    spot["input"] = input
    spot["output"] = output
    spot["stat"] = {}
    for spot_file in input["spot_files"]:
        with open(spot_file, "rU") as f:
            stat = f.readlines()[1].strip().split()
            stat_dict = {}
            stat_dict["total"] = float(stat[0])
            stat_dict["hotspot"] = float(stat[1])
            stat_dict["spot"] = float(stat[2])
            spot["stat"][spot_file] = stat_dict
    json_dump(spot)

def stat_peaks_number(input = {"peaks": ""}, output = {}):
    """ evaluate two passes hotspot and peaks merge regions number """
    data = open(input["peaks"]).readlines()
    json_dict = {"input": input, "output": output, "stat": len(data)}

def stat_reps(input={"peaks": ""},
             output={"json": ""}, param=None):
    ## overlap percentage between replicates
    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dump(result_dict)

def stat_cor(input={"correlation_R":"", "cor_pdf": "", "venn": "", },
             output={"json": ""}, param=None):
    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    with open(input["cor"], 'rU') as f:
        value = f.readline()
        values = value.strip().split()
    correlation_list = [float(i) for i in values[1:-1]]
    print(correlation_list)
    result_dict["stat"]["cor"] = correlation_list
    result_dict["stat"]["min_cor"] = min(correlation_list)
    result_dict["stat"]["judge"] = "Pass" if result_dict["stat"]["min_cor"] >= 0.6 else "Fail"
    result_dict["stat"]["cutoff"] = 0.6
    json_dump(result_dict)

def stat_promotor(input = {"ceas_r": ""}, output = {"json": ""}, param=None):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    promotor_percentage = []
    with open(input["ceas_r"], 'rU') as f:
        for line in f.readlines():
            if line.startswith("pie"):
                extract = re.findall("pie\(x=x,labels=c\(" + "(.*)" + "\)", line)
                content =  extract[0].split(",")
                promotor_percentage.append(content[0])
    json_dict["stat"]["genome"] = promotor_percentage[0]
    json_dict["stat"]["peaks"] = promotor_percentage[1]
    json_dump(json_dict)

def stat_conserv(input = "", output = "", param = None):
    json_dict = {"input": input, "output": output, "param": param, "stat": ""}
    with open(input, "rU") as f:
        score = f.read().strip()
        json_dict["stat"] = score
    json_dump(json_dict)

def stat_dhs(input={"pks_spot_bed": "", "dhs_peaks": ""}, output={"json": ""},
             param={}):
    """
    overlap with union DHS
    """
    peaks_info = {}
    peaks_info["total"] = len(open(input["pks_spot_bed"], 'r').readlines())
    peaks_info["dhs"] = len(open(input["dhs_peaks"], 'r').readlines())
    peaks_info['dhspercentage'] = peaks_info["dhs"] / peaks_info["total"]
    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    result_dict["stat"] = peaks_info
    result_dict["stat"]["cutoff"] = 0.8
    result_dict["stat"]["judge"] = "Pass" if result_dict["stat"]["dhspercentage"] >= 0.8 else "Fail"
    json_dump(result_dict)