#############################################################
# The first step for GCAP QC, add multiple format support
#    1. input fastq files, use fastqStatsAndSubsample
#    2. input BAM files, use sampleBam and fastqc to evaluate
#############################################################

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.sampling import sample_reads
from gcap.funcs.helpers import *

def _versatile_format(workflow, conf):
    """
    support SAM/BAM and BED reads files
    """
    for n, target in enumerate(conf.treatment_targets):
        if conf.seq_type.startswith("bam"):
            attach_back(workflow,
                ShellCommand(
                    "{tool} view -h {input[bam]} > {output[sam]}",
                    tool = "samtools",
                    input = {"bam": conf.treatment_bam[n]},
                    output = {"sam": target + "_all.sam"},
                    name = "samtools conversion"))
        elif conf.seq_type.startswith("sam"):
            attach_back(workflow,
                ShellCommand(
                    "{tool} -sf {input[sam]} {output[sam]}",
                    tool = "ln",
                    input = {"sam": conf.treatment_bam[n]}, ## actually sam files
                    output = {"sam": target + "_all.sam"}))
        elif conf.seq_type.startswith("bed"):
            attach_back(workflow,
                ShellCommand(
                    "{tool} -sf {input[bed]} {output[bed]}",
                    tool = "ln",
                    input = {"bed":conf.treatment_bed[n]}, ## actually bed files
                    output = {"bed":target + "_all.bed"}))

## parse fastqc result
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
    parse per sequence quality median from fastqc text file
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

## parse json files to latex reads length and sequence quality part
def seq_quality_doc(input = {"tex": "", "json": ""}, output = {}, param = {"reps": "", "se_samples": "", "pe_samples": ""}):
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
        param = {"section_name": "sequence_quality.py",
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

## fastqc command line
## fastqStatsAndSubsample command line
def seq_quality(workflow, conf, tex):
    """
    1a. use top 100k reads instead of sampled reads for fastqc evaluation
    at the same time, prepare input files
    """
    if conf.seq_type == "pe": ## PE
        sample_reads(workflow, conf, 100000, "fastq")
        for n, target in enumerate(conf.treatment_pair_data):
            for p in target:
                fastqc_run = attach_back(workflow,
                    ShellCommand(
                        "{tool} {input[fastq_sample]} --extract -t {param[threads]} -o {output[target_dir]}",
                        input= {"fastq_sample": p + "_100k.fastq"},
                        output={"target_dir": conf.target_dir,
                                "fastqc_summary": p + "_100k_fastqc/fastqc_data.txt"},
                        tool="fastqc",
                        param={"threads": 4}))
                fastqc_run.update(param=conf.items("fastqc"))
            attach_back(workflow, PythonCommand(stat_fastqc,
                input = {"fastqc_summaries": [ p + "_100k_fastqc/fastqc_data.txt" for p in target ]},
                output = {"json": conf.json_prefix + "_rep" + str(n+1) + "_fastqc.json"},
                param = {"samples": target }))
        attach_back(workflow, PythonCommand(
            seq_quality_doc,
            input = {"tex": tex, "json": [ conf.json_prefix + "_rep" + str(n+1) + "_fastqc.json" for n in range(len(conf.treatment_pair_data)) ]},
            output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
            param = {"seq_type": conf.seq_type, "reps": len(conf.treatment_pairs),
                     "pe_samples": [ target for target in conf.treatment_pair_data ]}))

    elif conf.seq_type == "se":
        sample_reads(workflow, conf, 100000, "fastq")
        for target in conf.treatment_targets:
            fastqc_run = attach_back(workflow,
                ShellCommand(
                    "{tool} {input[fastq_sample]} --extract -t {param[threads]} -o {output[target_dir]}",
                    input= {"fastq_sample": target +  "_100k.fastq"},
                    output={"target_dir": conf.target_dir,
                            "fastqc_summary": target + "_100k_fastqc/fastqc_data.txt"},
                    tool="fastqc",
                    param={"threads": 4},
                    name = "fastqc"))
            fastqc_run.update(param=conf.items("fastqc"))
        attach_back(workflow, PythonCommand(stat_fastqc,
            input = {"fastqc_summaries": [ t + "_100k_fastqc/fastqc_data.txt" for t in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_fastqc.json"},
            param = {"samples": conf.treatment_bases}))
        attach_back(workflow, PythonCommand(
            seq_quality_doc,
            input = {"tex": tex, "json": conf.json_prefix + "_fastqc.json"},
            output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
            param = {"seq_type": conf.seq_type, "reps": len(conf.treatment_pairs), "se_samples": conf.treatment_bases}))

    elif conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam"):
        _versatile_format(workflow, conf)
        sample_reads(workflow, conf, 100000, "sam")
        for target in conf.treatment_targets:
            attach_back(workflow,
                ShellCommand(
                    "{tool} view -bt {input[chrom_len]} {input[sam]} -o {param[tmp_bam]} && \
                    {tool} sort -m {param[max_mem]} {param[tmp_bam]} {param[output_prefix]}",
                    tool="samtools",
                    input={"sam": target + "_100k.sam", "chrom_len": conf.get_path("lib", "chrom_len")},
                    output={"bam": target + "_100k.bam"},
                    param={"tmp_bam": target + ".tmp.bam", "output_prefix": target + "_100k",
                           "max_mem": 5000000000},
                    name = "sampled sam to bam")) # Use 5G memory as default
        file_suffix = "_100k.bam"
        suffix = "_100k" ## output suffix
        ## for single end 100k fastq or PE, SE 100k sam
        for target in conf.treatment_targets:
            fastqc_run = attach_back(workflow,
                ShellCommand(
                    "{tool} {input[fastq_sample]} --extract -t {param[threads]} -o {output[target_dir]}",
                    input= {"fastq_sample": target +  file_suffix},
                    output={"target_dir": conf.target_dir,
                            "fastqc_summary": target + "%s_fastqc/fastqc_data.txt" % suffix},
                    tool="fastqc",
                    param={"threads": 4},
                    name = "fastqc"))
            fastqc_run.update(param=conf.items("fastqc"))
        ## get per sequence quality median
        attach_back(workflow, PythonCommand(stat_fastqc,
            input = {"fastqc_summaries": [ t + "%s_fastqc/fastqc_data.txt" % suffix for t in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_fastqc.json"},
            param = {"samples": conf.treatment_bases}))
        ## sequence quality latex load, reads length
        attach_back(workflow, PythonCommand(
            seq_quality_doc,
            input = {"tex": tex, "json": conf.json_prefix + "_fastqc.json"},
            output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
            param = {"seq_type": conf.seq_type, "reps": len(conf.treatment_pairs), "se_samples": conf.treatment_bases})) ## single bam samples

    elif conf.seq_type.startswith("bed"): ## not evaluate sequence quality
        _versatile_format(workflow, conf)