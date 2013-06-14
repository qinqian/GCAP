from configparser import ConfigParser
import json
import os
import random
import re
import sys

from jinja2 import Environment, FileSystemLoader
from samflow.command import AbstractCommand
from samflow.workflow import attach_back


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


def count_in_million(x):
    if type(x) == int:
        return str(round(x/1000000, 1)) + "M"

def decimal_to_latex_percent(dec):
    if type(dec) == float:
        return str(round(dec*100, 1))    + "\%"

def underline_to_space(x):
    if type(x) == str:
        return x.replace("_", " ")
    return x

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


def begin_doc(input = "", output = {"begin": "", "header": ""}, param = {"rep": "", "id": ""}):
    ## begin document
    begin_latex = JinjaTemplateCommand(
        template = input,
        name = "begin",
        param = {"section_name": "begin",
                 "render_dump": output["begin"],
                 "dataset_id": param["id"]})
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

def DHS_doc(input = {"json": "", "tex": ""}, output = {"latex": ""}, param = {"reps": "", "samples": ""}):
    data = json_load(input["json"])["stat"]
    dhs = []

    for s in param["samples"]:
        dhs.append(decimal_to_latex_percent(data[s]["dhspercentage"]))

    DHS_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "DHS",
        param = {"section_name": "DHS",
                 "render_dump": output["latex"],
                 "DHS": dhs,
                 "reps": param["reps"]})
    template_dump(DHS_latex)

def reads_doc(input = {"tex": "", #"json": "",
                       "json_autosome": ""},
              output = {"raw": "", "mapping": ""}, param = {"seq_type": "", "reps": "", "samples": ""}):
    """
    these values will be replaced by Jim's codes,
    currently by bowtie and fastqc,
    use autosome reads mappable ratio
    """
    json_au = json_load(input["json_autosome"])["stat"]  ## autosome
    ## json = json_load(input["json"])["stat"]
    ## use correct order

    if param["seq_type"] == "se":
        raw =  [ count_in_million(json_au[s]["total"]) for s in param["samples"] ]
        total = sum([int(json_au[j]["total"]) for j in json_au])
        mapped_number = [ int(json_au[s]["map"]) for s in param["samples"] ]
    else:
        raw =  [ "%s, %s" % (count_in_million(int(json_au[s]["total"])/2), count_in_million(int(json_au[s]["total"])/2)) for s in param["samples"] ]
        total = sum([int(json_au[j]["total"]) for j in json_au])
        mapped_number = [int(json_au[s]["map"]) for s in param["samples"]]

    ## use autosome ratio
    mapped_rate = [ decimal_to_latex_percent(json_au[s]["ratio"]) for s in param["samples"] ]

    ## raw reads doc
    reads_latex = JinjaTemplateCommand(
        template = input['tex'],
        name = "raw_reads",
        param = {"section_name": "raw_reads",
                 "render_dump": output["raw"],
                 "reads": raw,
                 "reps": param["reps"],
                 "combo": count_in_million(total)})
    template_dump(reads_latex)

    ## mapping status doc
    mapping_latex = JinjaTemplateCommand(
        template = input['tex'],
        name = "reads mapping",
        param = {"section_name": "mapping",
                 "render_dump": output["mapping"],
                 "map": mapped_rate, "combo": count_in_million(sum(mapped_number)),
                 "reps": param["reps"]})

    template_dump(mapping_latex)

def fastqc_doc(input = {"tex": "", "json": ""}, output = {}, param = {"reps": "", "se_samples": "", "pe_samples": ""}):
    if type(input["json"]) == str:
        json = json_load(input['json'])
        seq = [ json["stat"][s]["median"] for s in param["se_samples"] ]
        len = [ str(json['stat'][s]['sequence_length']) + " " + param["seq_type"].upper() for s in param["se_samples"] ]
    elif type(input["json"]) == list:
        seq = []
        len = []

        for j, s in zip(input["json"], param["pe_samples"]):
            pair_seq = []
            pair_len = []
            data = json_load(j)

            for d in s:
                pair_seq.append(str(data['stat'][d]['median']))
                pair_len.append(str(data['stat'][d]['sequence_length']) + " " + param["seq_type"].upper())

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

def redundancy_doc(input = {"tex": "", "json": ""}, output = {"redun": ""}, param = {"reps": "", "samples": ""}):

    redun = json_load(input["json"])["stat"]
    data = []

    for s in param["samples"]:
        data.append(decimal_to_latex_percent(float(redun[s])))

    print(data)
    redun_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "Library complexity",
        param = {"section_name": "redundancy",
                 "render_dump": output["redun"],
                 "redun": data,
                 "reps": param["reps"]})

    template_dump(redun_latex)

def stat_redun(input = {"picard": "", "built_count": ""}, output = {"json": ""}, param = {"samples": "", "format": ""}):
    """
    this will be replaced by Gifford's codes,
    currently picard
    """
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    if param["format"] != "bed":
        for f, s in zip(input["picard"], param["samples"]):
            data = open(f).readlines()[7]
            d = data.split("\t")[7]
            json_dict["stat"][s] = d
    else:
        for b, s in zip(input["built_count"], param["samples"]):
            total_loc = open(b[0]).read().strip().split()[1]
            uniq_loc = open(b[1]).read().strip().split()[1]
            redun_ratio = 1 - float(uniq_loc) / total_loc
            json_dict["stat"][s] = redun_ratio
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
                frag_size = re.findall("c\((.*\d+.+)\)", line)
                values["frag"] = frag_size[0].split(",")
    return values

def stat_frag_std(input = {"r": "", "insert": ""}, output = {"json": "", "r": ""}, param = {"samples": "", "frag_tool": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    if param["frag_tool"] == "macs2":
        ## single end
        for rin, rout, s in zip(input["r"], output["r"], param["samples"]):
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
            json_dict["stat"][s] = "mean %s, sd %s" % (max(map(int, values["frag"])), int(float(std)))

    elif param["frag_tool"] == "picard":
        for i, s in zip(input["insert"], param["samples"]):
            with open(i, "rU") as dup:
                data = dup.readlines()[7]
                # print(data.split("\t")[4])
                # print(data.split("\t")[5])
            data = data.strip().split("\t")
            json_dict["stat"][s] = str(int(float(data[4]))) + ", sd " + str(int(float(data[5])))
    json_dump(json_dict)

def frag_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"reps": "", "samples": ""}):

    frag = json_load(input["json"])["stat"]
    data = []
    for s in param["samples"]:
        data.append(frag[s])

    frag_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "fragment",
        param = {"section_name": "fragment",
                 "render_dump": output["latex"],
                 "frag": data,
                 "reps": param["reps"]})
    template_dump(frag_latex)

def stat_peaks(input = {"peaks": {"all_peaks": "", "5M_spot": "", "combo_peaks": ""}}, output = {"json"}, param = {"tool": "", "samples": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}

    json_dict["stat"]["all_peaks"] = {}

    # d, narrow peaks
    for a, s in zip(input["peaks"]["all_peaks"], param["samples"]):
        with open(a) as f:
           json_dict["stat"]["all_peaks"][s] = len(f.readlines())

    json_dict["stat"]["5M_spot"] = {}

    # b, hotspot minimally threshold
    for spot, sam in zip(input["peaks"]["5M_spot"], param["samples"]):
        with open(spot) as f:
            json_dict["stat"]["5M_spot"][sam] = len(f.readlines()[1:])

    json_dict["stat"]["combo"] = {}

    if len(input["peaks"]["all_peaks"]) >= 2:
        with open(input["peaks"]["combo"]) as combo:
            json_dict["stat"]["combo"]= len(combo.readlines())
    print(json_dict)

    json_dump(json_dict)

def peaks_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"samples": "", "reps": ""}):
    """ use d peaks( narrow peaks ) for all reads evaluation
     use b Hotspot for 5M reads evaluation and replicates consistency(correlation and overlap)
     use 5M narrow peaks for DHS evaluation
    """
    data = json_load(input["json"])
    peaks_5M = []
    peaks_all = []

    for s in param["samples"]:
        peaks_5M.append(data["stat"]["5M_spot"][s])

    for a in param["samples"]:
        peaks_all.append(data["stat"]["all_peaks"][a])

    if data["stat"]["combo"]:
        combo = data["stat"]["combo"]
        peaks_latex = JinjaTemplateCommand(
            template = input["tex"],
            name = "peaks_calling",
            param = {"section_name": "peaks_calling",
                     "render_dump": output["latex"],
                     "spot_5M": peaks_5M,
                     "peaks_all":peaks_all, "combo": combo,
                     "tool": data["param"]["tool"],
                     "reps": param["reps"]})

    else:
        peaks_latex = JinjaTemplateCommand(
            template = input["tex"],
            name = "peaks_calling",
            param = {"section_name": "peaks_calling",
                     "render_dump": output["latex"],
                     "spot_5M": peaks_5M,
                     "peaks_all": peaks_all,
                     "tool": data["param"]["tool"],
                     "reps": param["reps"]})
    template_dump(peaks_latex)

def promotor_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"samples": "", "reps": ""}):

    stat = json_load(input["json"])["stat"]
    genome = stat["genome_promotor_percentage"]


    peaks = [ stat["promotor_percentage"][s] for s in param["samples"] ]

    promotor_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "promotor",
        param = {"section_name": "promotor",
                 "render_dump": output["latex"],
                 "promotor": peaks,
                 "genome": genome,
                 "reps": param["reps"]})

    template_dump(promotor_latex)


def single_end_fastq_sampling(input = {"fastq": ""}, output = {"fastq_sample": ""}, param = {"random_number": ""}):
    num_lines = sum(1 for _ in open(input["fastq"])) / 4

    rand_nums = sorted([random.randint(0, num_lines - 1) for _ in range(param["random_number"])])

    fastq = open(input["fastq"],"rU")
    fastq_sample = open(output["fastq_sample"], "w")

    cur_num = - 1
    written = 0

    for rand_num in rand_nums:
        while cur_num < rand_num:
            for i in range(4):
                fastq.readline()
            cur_num += 1

        for i in range(4):
            fastq_sample.write(fastq.readline())
        cur_num += 1
        written += 1
    assert written == param["random_number"]

    fastq_sample.close()
    fastq.close()

def pair_end_fastq_sampling(input = {"fastq": ""}, output = {"fastq_sample": ""}, param = {"random_number": ""}):
    def write_random_records(fqa, fqb, suba, subb, N=100000):
        records = sum(1 for _ in open(fqa)) / 4
        rand_records = sorted([random.randint(0, records - 1) for _ in range(N)])

        fha, fhb = open(fqa),  open(fqb)
        suba, subb = open(suba, "w"), open(subb, "w")
        rec_no = - 1
        written = 0

        for rr in rand_records:
            while rec_no < rr:
                rec_no += 1
                for i in range(4): fha.readline()
                for i in range(4): fhb.readline()
            for i in range(4):
                suba.write(fha.readline())
                subb.write(fhb.readline())
            rec_no += 1
            written += 1
        assert written == N
        suba.close()
        subb.close()
        fha.close()
        fhb.close()
    write_random_records(input["fastq"][0], input["fastq"][1], output["fastq_sample"][0], output["fastq_sample"][1], param["random_number"])

def sampling_sam(input = {"sam": ""}, output = {"sam_sample": ""}, param = {"random_number": "", "map_or_unmap": "both", "se_or_pe": ""}):
    """
    sampling SAM or BED reads files,
    need to add map_or_unmap
    """
    num_lines = sum(1 for _ in open(input["sam"]))
    header_num = 0
    header = []
    with open(input["sam"]) as f:
        for line in f:
            if line.startswith("@"):
                header_num += 1
                header.append(line)
            else:
                break
    print(header_num)

    if param["map_or_unmap"] == "both":
        if param["se_or_pe"] == "se":
            rand_nums = sorted([random.randint(header_num, num_lines - 1) for _ in range(param["random_number"])])
            print(len(rand_nums))

            cur_num = -1
            written = 0

            with open(output["sam_sample"], "w") as fout:
                with open(input["sam"], "rU") as fin:
                    for rand_num in rand_nums:
                        while cur_num < rand_num:
                            cur_num+=1
                            data = fin.readline()
                            if data.startswith("@"):
                                fout.write(data)

                        fout.write(fin.readline())
                        cur_num += 1
                        written += 1
            assert  written == param["random_number"]
        elif param["se_or_pe"] == "pe":
            random_pe = []
            cur_num = -1
            written = 0

            ## large memory version 1, faster
            for _ in range(int(param["random_number"]/2)):
                num = random.choice(range(header_num, num_lines, 2))
                random_pe.append(num)
                random_pe.append(num+1)
            rand_nums = random_pe
            data = open(input["sam"]).readlines()

            with open(output["sam_sample"], 'w') as f:
               f.write("".join(header))
               for i in rand_nums:
                   f.write(data[i])

            ## small memory version 2, slower and simpler
            ## most of pair end reads paired, need to modify
#            for _ in range(int(param["random_number"]/2)):
#               num = random.choice(range(header_num, num_lines, 2))
#               random_pe.append(num)
#            rand_nums = sorted(random_pe)
#            with open(output["sam_sample"], "w") as fout:
#                with open(input["sam"], "rU") as fin:
#                    for rand_num in rand_nums:
#                        while cur_num < rand_num:
#                            cur_num+=1
#                            data = fin.readline()
#                            if data.startswith("@"):
#                                fout.write(data)
#                        fout.write(fin.readline())
#                        fout.write(fin.readline())
#                        cur_num += 1
#                        written += 2
#                assert  written == param["random_number"]
    else: pass

def autosome_map(input = {"count": ""}, output = {"json": ""}, param = {"samples": ""}):
    """
    use autosome mapping ratio to replace all reads mapping ratio
    """
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    for i, s in zip(input["count"], param["samples"]):
        json_dict["stat"][s] = {}
        with open(i, 'rU') as f:
            data = f.read().strip().split("\t")
            no_chrM_tag = data[0]
            total = data[1]

        json_dict["stat"][s]["ratio"] = round(float(no_chrM_tag) / float(total), 3)
        json_dict["stat"][s]["map"] = no_chrM_tag
        json_dict["stat"][s]["total"] = total

    json_dump(json_dict)

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
                param= {"samples": ""}):
    """
    per sequence quality
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    stat = {}
    for a_summary, a_id in zip(input["fastqc_summaries"], param["samples"]):
        print(a_summary, a_id)
        parsed = fastqc_parse(input=a_summary)
        stat[a_id] = {}
        stat[a_id]["median"] = parsed["median"]
        stat[a_id]["cutoff"] = 25
        stat[a_id]["sequence_length"] = parsed["sequence_length"]

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
              param = {"K": "", "FDRS": "0.01", "species": "", "keepdup": "T"}):
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
    setting(cf, "_DUPOK_", param["keepdup"])         ## keep duplicates
    setting(cf, "_CHKCHR_", "chrX")     ## To check if data contains this chromosome
    setting(cf, "_CHECK_", "F")         ## change to F in case that hotspot error for no wavelets peaks, other method
    cf.write(open(output["conf"], "w"))

def stat_spot_on_replicates(input = {"spot_files": ""}, output = {"json": ""},
                            param = {"samples": ""}):
    spot = {}
    spot["input"] = input
    spot["output"] = output
    spot["stat"] = {}
    for spot_file, s in zip(input["spot_files"], param["samples"]):
        with open(spot_file, "rU") as f:
            stat = f.readlines()[1].strip().split()
            stat_dict = {}
            stat_dict["total"] = float(stat[0])
            stat_dict["hotspot"] = float(stat[1])
            stat_dict["spot"] = float(stat[2])
            spot["stat"][s] = stat_dict
    json_dump(spot)

def spot_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"samples": "", "reps": ""}):
    data = json_load(input["json"])["stat"]
    spot = []
    for s in param["samples"]:
        spot.append(decimal_to_latex_percent(data[s]["spot"]))

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
        json_dict["stat"]["overlap"][rep[1]] = "%s" % (decimal_to_latex_percent(float(rep_overlap) / union_num))

    json_dict["stat"]["cor"] = {}
    for rep in input["5M_cor"]:
        with open(rep[0]) as f:
            rep_cor = f.read().strip().split()[0]

        json_dict["stat"]["cor"][rep[1]] = "%s" % rep_cor

    json_dump(json_dict)


def reps_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {}):
    data = json_load(input["json"])
    cor = data['stat']["cor"]
    cor_d = []
    for i in cor:
        cor_d.append(str(round(float(cor[i]), 2)))

    overlap = data["stat"]["overlap"]
    overlap_d = []
    for i in overlap:
        overlap_d.append(overlap[i])


    rep_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "replicates",
        param = {"section_name": "replicates",
                 "render_dump": output["latex"],
                 "overlap": overlap_d, "cor": cor_d})

    template_dump(rep_latex)

def stat_promotor(input = {"peaks_promotor": "", "peaks": "", "promotor": "", "mappable": ""}, output = {"json": ""}, param={"samples": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}

    json_dict["stat"]["promotor_percentage"] = {}
    for overlap, peaks, s in zip(input["peaks_promotor"], input["peaks"], param["samples"]):
        with open(overlap) as p:
            overlap_num = p.read().strip().split()[0]
        with open(peaks) as f:
            peaks_num = f.read().strip().split()[0]

        json_dict["stat"]["promotor_percentage"][s] = decimal_to_latex_percent(float(overlap_num) / float(peaks_num))

    with open(input["promotor"]) as g:
        promotor = g.read().strip().split()[0]

    with open(input["mappable"]) as m:
        mappable = len(m.readlines())

    json_dict["stat"]["genome_promotor_percentage"] = decimal_to_latex_percent(float(promotor) / mappable)

    json_dump(json_dict)


def stat_conserv(input = {"phastcon": ""}, output = {"json": ""}, param = {"sample": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}

    for f, s in zip(input["phastcon"], param["sample"]):
        print(f, s)
        with open(f, "rU") as phast:
            score = phast.read().strip()
            json_dict["stat"][s] = score

    json_dump(json_dict)

def conserv_doc(input = {"json": "", "tex": ""}, output = {"latex": ""}, param = {"reps": "", "sample": ""}):
    data = json_load(input["json"])["stat"]
    conserv = []
    for s in param["sample"]:
        conserv.append(round(float(data[s]), 2))

    conservation_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "conservation",
        param = {"section_name": "conservation",
                 "render_dump": output["latex"],
                 "conservation": conserv,
                 "reps": param["reps"]})
    template_dump(conservation_latex)

def stat_dhs(input={"pks_spot_bed": "", "dhs_peaks": ""}, output={"json": ""},
             param={"samples": ""}):
    """
    overlap with union DHS
    """
    peaks_info = {}
    for b, d, s in zip(input["pks_spot_bed"], input["dhs_peaks"], param["samples"]):
        print(b, d, s)

        peaks_info[s] = {}
        peaks_info[s]["total"] = len(open(b, 'r').readlines())
        peaks_info[s]["dhs"] = len(open(d, 'r').readlines())
        peaks_info[s]['dhspercentage'] = peaks_info[s]["dhs"] / peaks_info[s]["total"]
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict["stat"] = peaks_info
    json_dump(json_dict)
