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


def begin_doc(input = "", output = {"begin": "", "header": ""}, param = {"rep": ""}):
    ## begin document
    begin_latex = JinjaTemplateCommand(
        template = input,
        name = "begin",
        param = {"section_name": "begin",
                 "render_dump": output["begin"]})
    template_dump(begin_latex)

    header_latex = JinjaTemplateCommand(
        template = input,
        name = "header",
        param = {"section_name": "header",
                 "render_dump": output["header"],
                 "reps": param["rep"]})
    template_dump(header_latex)

def end_doc(input = {"tex": ""}, output = {"table_end": "", "doc_end": ""}, param = {}):
    """
    table end and latex end
    """
    table_end_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "end",
        param = {"section_name": "table_end",
                 "render_dump": output["table_end"]})

    template_dump(table_end_latex)

    end_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "end",
        param = {"section_name": "end",
                 "render_dump": output["doc_end"]})
    template_dump(end_latex)

def DHS_doc(input = {"json": "", "tex": ""}, output = {"latex": ""}, param = {"reps": ""}):
    data = json_load(input["json"])["stat"]
    print(data)
    dhs = []
    for s in data:
        dhs.append(decimal_to_latex_percent(data[s]["dhspercentage"]))

    DHS_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "DHS",
        param = {"section_name": "DHS",
                 "render_dump": output["latex"],
                 "DHS": dhs,
                 "reps": param["reps"]})
    template_dump(DHS_latex)

def reads_doc(input = {"tex": "", "json": ""}, output = {"raw": "", "mapping": ""}, param = {"seq_type": "", "reps": ""}):
    """
    these values will be replaced by Jim's codes,
    currently by bowtie and fastqc
    """
    json = json_load(input["json"])["stat"]

    if param["seq_type"] == "se":
        raw =  [ json[j]["total_reads"] for j in json ]
        mapped_number = [ json[j]["mappable_reads"] for j in json ]
    else:
        raw =  [ json[j]["total_reads"] * 2 for j in json ]
        mapped_number = [ json[j]["mappable_reads"] * 2 for j in json ]

    mapped_rate = [ json[j]["mappable_rate"] for j in json ]

    print(param["reps"])

    reads_latex = JinjaTemplateCommand(
        template = input['tex'],
        name = "raw_reads",
        param = {"section_name": "raw_reads",
                 "render_dump": output["raw"],
                 "reads": raw,
                 "reps": param["reps"],
                 "combo": sum(raw)})

    template_dump(reads_latex)

    mapping_latex = JinjaTemplateCommand(
        template = input['tex'],
        name = "reads mapping",
        param = {"section_name": "mapping",
                 "render_dump": output["mapping"],
                 "map": mapped_rate, "combo": sum(mapped_number),
                 "reps": param["reps"]})

    template_dump(mapping_latex)

def fastqc_doc(input = {"tex": "", "json": ""}, output = {}, param = {"reps": ""}):
    if type(input["json"]) == str:
        json = json_load(input['json'])
        seq = [ json["stat"][d]["median"] for d in json["stat"] ]
        len = [ str(json['stat'][d]['sequence_length']) + param["seq_type"].upper() for d in json['stat'] ]
    elif type(input["json"]) == list:
        seq = []
        len = []
        for j in input["json"]:
            pair_seq = []
            pair_len = []
            data = json_load(j)
            for d in data["stat"]:
                pair_seq.append(str(data['stat'][d]['median']))
                pair_len.append(str(data['stat'][d]['sequence_length']) + param["seq_type"].upper())

            seq.append(','.join(pair_seq))
            len.append(','.join(pair_len))

    seq_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "sequence quality",
        param = {"section_name": "sequence_quality",
                 "render_dump": output["seq"],
                 "seq_quality": seq,
                 "reps": param["reps"]})

    len_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "reads_len",
        param = {"section_name": "reads_length",
                 "render_dump": output["len"],
                 "reads_len": len,
                 "reps": param["reps"]})

    template_dump(seq_latex)
    template_dump(len_latex)

def redundancy_doc(input = {"tex": "", "json": ""}, output = {"redun": ""}, param = {"reps": ""}):

    redun = json_load(input["json"])["stat"]
    data = []
    for d in redun:
        print(d)
        data.append(redun[d])

    redun_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "Library complexity",
        param = {"section_name": "redundancy",
                 "render_dump": output["redun"],
                 "redun": data,
                 "reps": param["reps"]})

    template_dump(redun_latex)

def stat_redun(input = {"picard": ""}, output = {"json": ""}, param = {}):
    """
    this will be replaced by Gifford's codes,
    currently picard
    """
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    for f in input["picard"]:
        data = open(f).readlines()[7]
        d = data.split("\t")[7]
        json_dict["stat"][f] = d

    json_dump(json_dict)

def get_size(rscript):
    values = {}
    with open(rscript) as model:
        for line in model:
            if line.startswith("p <- "):
                values["positive"] = line
            if line.startswith("m <- "):
                values["minus"] = line
            if line.startswith("x <- "):
                values["x"] = line
            if line.startswith("altd"):
                frag_size = re.findall("c\((\d+)\)", line)
                values["frag"] = frag_size
    return values

def stat_frag_std(input = {"r": ""}, output = {"json": "", "r": ""}, param = {}):
    json_dict = {"input": input, "output": output, "param": param, "stat": []}
    for rin, rout in zip(input["r"], output["r"]):
        print(rin)
        print(rout)
        values = get_size(rin)
        with open(rout, 'w') as f:
            f.write(values['positive'])
            f.write(values['minus'])
            f.write(values['x'])
            f.write("p.expect = sum(x * p/100) \n")
            f.write("m.expect = sum(x * m/100) \n")
            f.write("p.sd = sqrt(sum(((x-p.expect)^2)*p/100)) \n")
            f.write("m.sd = sqrt(sum(((x-m.expect)^2)*m/100)) \n")
            f.write("cat(paste((p.sd + m.sd)/2, '\n')) \n")

        std = os.popen("Rscript %s" % rout).read().strip()
        json_dict["stat"].append("mean %s, sd %s" % (std, values["frag"][0]))
    json_dump(json_dict)

def frag_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"reps": ""}):

    frag = json_load(input["json"])["stat"]

    frag_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "fragment",
        param = {"section_name": "fragment",
                 "render_dump": output["latex"],
                 "frag": frag,
                 "reps": param["reps"]})
    template_dump(frag_latex)

def stat_peaks(input = {"peaks": {"all": "", "5M": "", "combo": ""}}, output = {"json"}, param = {"tool": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}

    json_dict["stat"]["all"] = {}

    for a in input["peaks"]["all"]:
        with open(a) as f:
           json_dict["stat"]["all"][a] = len(f.readlines())

    json_dict["stat"]["5M"] = {}
    for s in input["peaks"]["5M"]:
        with open(a) as f:
            json_dict["stat"]["5M"][s] = len(f.readlines())

    json_dict["stat"]["combo"] = {}
    if len(input["peaks"]["all"]) >= 2:
        with open(input["peaks"]["combo"]) as combo:
            json_dict["stat"]["combo"]= len(combo.readlines())
    json_dump(json_dict)

def peaks_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"samples": "", "reps": ""}):
    data = json_load(input["json"])
    peaks_5M = data["stat"]["5M"].values()
    peaks_all = data["stat"]["all"].values()

    if "combo" in data["stat"]:

        combo = data["stat"]["combo"]
        peaks_latex = JinjaTemplateCommand(
            template = input["tex"],
            name = "peaks_calling",
            param = {"section_name": "peaks_calling",
                     "render_dump": output["latex"],
                     "peaks_5M": peaks_5M,
                     "peaks_all":peaks_all, "combo": combo,
                     "tool": data["param"]["tool"],
                     "reps": param["reps"]})

    else:
        peaks_latex = JinjaTemplateCommand(
            template = input["tex"],
            name = "peaks_calling",
            param = {"section_name": "peaks_calling",
                     "render_dump": output["latex"],
                     "peaks_5M": peaks_5M,
                     "peaks_all": peaks_all,
                     "tool": data["param"]["tool"]})
    template_dump(peaks_latex)



def promotor_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"samples": "", "reps": ""}):

    stat = json_load(input["json"])["stat"]
    genome = stat["genome_promotor_percentage"]
    peaks = stat["promotor_percentage"].values()

    add_genome = lambda x: "peaks: " + x + ", genome: " + genome
    peaks_genome = map(add_genome, peaks)

    promotor_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "promotor",
        param = {"section_name": "promotor",
                 "render_dump": output["latex"],
                 "promotor": peaks_genome,
                 "reps": param["reps"]})

    template_dump(promotor_latex)


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
                       output = {"json": ""}, param = {"samples": "", "species": "", "id": "", "correct_species": ""}):
    ## if library is contaminated, break
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

    for s in library_contamination["value"]:
        # s for different samples
        for i in param["species"]:
            # i for species, check each species
            if i != param["correct_species"]:
                if library_contamination["value"][s][param["correct_species"]] > library_contamination["value"][s][i] :
                    print("Library is correct !")
                else:
                    print("Library may be contaminated !")
                    sys.exit(1)

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

def stat_spot_on_replicates(input = {"spot_files": ""}, output = {"json": ""},
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

def spot_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"sample": "", "reps": ""}):
    data = json_load(input["json"])["stat"]
    spot = []
    for d in data:
        spot.append(data[d]["spot"])

    spot_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "peaks_calling",
        param = {"section_name": "spot",
                 "render_dump": output["latex"],
                 "spot": spot,
                 "reps": param["reps"]})

    template_dump(spot_latex)

def stat_reps(input={"5M_overlap": "", "5M_cor": "", "union": ""},
             output={"json": ""}, param=None):
    ## overlap percentage between replicates
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    with open(input["union"]) as f:
        union_num = len(f.readlines())

    json_dict["stat"]["overlap"] = {}
    for rep in input["5M_overlap"]:
        rep_overlap = len(open(rep[0]).readlines())
        json_dict["stat"]["overlap"][rep[1]] = "%s %s, %s" % (rep[1], rep_overlap, decimal_to_latex_percent(float(rep_overlap) / union_num))

    json_dict["stat"]["cor"] = {}
    for rep in input["5M_cor"]:
        with open(rep[0]) as f:
            rep_cor = f.read().strip().split()[2]
        json_dict["stat"]["cor"][rep[1]] = "%s %s" % (rep[1], rep_cor)

    json_dump(json_dict)


def reps_doc(input = {"tex": "", "json": ""}, output = {}, param = {}):
    data = json_load(input["json"])
    cor = data['stat']["cor"].values()
    overlap = data["stat"]["overlap"].values()
    rep_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "replicates",
        param = {"section_name": "replicates",
                 "render_dump": "rep.tex",
                 "overlap": overlap, "cor": cor})
    template_dump(rep_latex)

def stat_promotor(input = {"peaks_promotor": "", "peaks": "", "promotor": "", "mappable": ""}, output = {"json": ""}, param=None):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}

    json_dict["stat"]["promotor_percentage"] = {}
    for overlap, peaks in zip(input["peaks_promotor"], input["peaks"]):
        with open(overlap) as p:
            overlap_num = p.read().strip().split()[0]
        with open(peaks) as f:
            peaks_num = f.read().strip().split()[0]

        json_dict["stat"]["promotor_percentage"][peaks] = decimal_to_latex_percent(float(overlap_num) / float(peaks_num))

    with open(input["promotor"]) as g:
        promotor = g.read().strip().split()[0]

    with open(input["mappable"]) as m:
        mappable = len(m.readlines())

    json_dict["stat"]["genome_promotor_percentage"] = decimal_to_latex_percent(float(promotor) / mappable)

    json_dump(json_dict)


def stat_conserv(input = {"phastcon": ""}, output = {"json": ""}, param = {"sample": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    print(input["phastcon"])
    print(param["sample"])
    for f, s in zip(input["phastcon"], param["sample"]):
        with open(f, "rU") as phast:
            score = phast.read().strip()

            json_dict["stat"][s] = score

    json_dump(json_dict)


def conserv_doc(input = {"json": "", "tex": ""}, output = {"latex": ""}, param = {"reps": ""}):
    data = json_load(input["json"])["stat"]
    conserv = data.values()

    conservation_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "conservation",
        param = {"section_name": "conservation",
                 "render_dump": output["latex"],
                 "conservation": conserv,
                 "reps": param["reps"]})
    template_dump(conservation_latex)

def stat_dhs(input={"pks_spot_bed": "", "dhs_peaks": ""}, output={"json": ""},
             param={}):
    """
    overlap with union DHS
    """
    peaks_info = {}
    for b, d in zip(input["pks_spot_bed"], input["dhs_peaks"]):
        peaks_info[b] = {}
        peaks_info[b]["total"] = len(open(b, 'r').readlines())
        peaks_info[b]["dhs"] = len(open(d, 'r').readlines())
        peaks_info[b]['dhspercentage'] = peaks_info[b]["dhs"] / peaks_info[b]["total"]

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict["stat"] = peaks_info
    json_dump(json_dict)
