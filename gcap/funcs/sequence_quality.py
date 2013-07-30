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
        ## link to uniform interface
        if conf.seq_type.startswith("bam"):
            attach_back(workflow,
                ShellCommand(
                    "{tool} -sf {input[bam]} {output[bam]} && \
                    samtools view -h {input[bam]} > {output[sam]}",
                    tool = "ln",
                    input = {"bam": conf.treatment_bam[n]},
                    output = {"bam": target + "_all.bam",
                              "sam": target + "_all.sam"}, ## for pipeline-scripts/autosome* for mapping statistics
                    name = "link"))
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

            if conf.peakcalltool == "hotspot": ## prepare bed.starch for hotspot peaks calling
                attach_back(workflow,
                    ShellCommand(
                        "sort-bed {input[bed]} | starch - > {output[starch]}",
                        tool = "ln",
                        input = {"bed":conf.treatment_bed[n]}, ## actually bed files
                        output = {"starch":conf.hotspot_starch_input[n] + "_all.bed.starch"}))

## parse fastqStatsAndSubsample result
def stat_fastqStat(input = {"seq": ""}, output = {"json": ""}, param = {"samples": "", "seq_type": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    if param["seq_type"] == "pe":
        for i, s in zip(input['seq'], param['samples']):
            with open(i[0]) as f:
                data = f.readlines()
                json_dict["stat"][s + "_pair1"] = {}
                json_dict["stat"][s + "_pair1"]['count'] = data[0].strip().split()[1]
                json_dict["stat"][s + "_pair1"]['quality'] = round(float(data[6].strip().split()[1]), 1)
                json_dict["stat"][s + "_pair1"]['std'] = data[7].strip().split()[1]
                json_dict["stat"][s + "_pair1"]['len'] = data[2].strip().split()[1]
            with open(i[1]) as f:
                data = f.readlines()
                json_dict["stat"][s + "_pair2"] = {}
                json_dict["stat"][s + "_pair2"]['count'] = data[0].strip().split()[1]
                json_dict["stat"][s + "_pair2"]['quality'] = round(float(data[6].strip().split()[1]), 1)
                json_dict["stat"][s + "_pair2"]['std'] = data[7].strip().split()[1]
                json_dict["stat"][s + "_pair2"]['len'] = data[2].strip().split()[1]
    if param["seq_type"] in ["se", "bam", "sam"]:
        for i, s in zip(input['seq'], param['samples']):
            with open(i) as f:
                data = f.readlines()
                json_dict["stat"][s] = {}
                json_dict["stat"][s]['count'] = data[0].strip().split()[1]
                json_dict["stat"][s]['quality'] = data[6].strip().split()[1]
                json_dict["stat"][s]['std'] = data[7].strip().split()[1]
                json_dict["stat"][s]['len'] = data[2].strip().split()[1]
    json_dump(json_dict)

## parse json files to latex reads length and sequence quality part
def seq_quality_doc(input = {"tex": "", "json": ""}, output = {}, param = {"reps": "", "se_samples": "", "pe_samples": "", "seq_type": ""}):
    data = json_load(input["json"])['stat']
    seq = []
    len = []
    if param["seq_type"] == "pe":
        for s in param["pe_samples"]:
            seq.append(str(data[s + "_pair1"]['quality']) + ',' + str(data[s + "_pair2"]['quality']))
            len.append(data[s + "_pair1"]['len'] + ',' + data[s + "_pair2"]['len'])
    else:
        for s in param["se_samples"]:
            seq.append(data[s]['quality'])
            len.append(data[s]['len'])

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

## fastqStatsAndSubsample command line
def seq_quality(workflow, conf, tex):
    """
    1a. use top 100k reads instead of sampled reads for fastqc evaluation
    at the same time, prepare input files
    """
    if conf.seq_type == "pe": ## PE
        sample_reads(workflow, conf, 100000, "fastq")
        attach_back(workflow, PythonCommand(stat_fastqStat,
            input = {"seq": [ [ p + "_100k.seq" for p in target ] for target in conf.treatment_pair_data ]},
            output = {"json": conf.json_prefix + "_seq_quality.json"},
            param = {"samples": conf.treatment_bases, "seq_type": conf.seq_type}))
        attach_back(workflow, PythonCommand(
            seq_quality_doc,
            input = {"tex": tex, "json": conf.json_prefix + "_seq_quality.json"},
            output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
            param = {"seq_type": conf.seq_type, "reps": len(conf.treatment_pairs),
                     "pe_samples": conf.treatment_bases}))

    elif conf.seq_type == "se":
        sample_reads(workflow, conf, 100000, "fastq")
        attach_back(workflow, PythonCommand(stat_fastqStat,
            input = {"seq": [ t + "_100k.seq" for t in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_seq_quality.json"},
            param = {"samples": conf.treatment_bases,
                     "seq_type": conf.seq_type}))
        attach_back(workflow, PythonCommand(
            seq_quality_doc,
            input = {"tex": tex, "json": conf.json_prefix + "_seq_quality.json"},
            output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
            param = {"seq_type": conf.seq_type, "reps": len(conf.treatment_pairs), "se_samples": conf.treatment_bases}))

    elif conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam"):
        _versatile_format(workflow, conf)
        if conf.seq_type.startswith("bam"):
            sample_reads(workflow, conf, 110000, "bam")
        else:
            sample_reads(workflow, conf, 110000, "sam")
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
            attach_back(workflow,
                ShellCommand(
                    "{tool} -i {input[bam]} -fq {output[fastq]}",
                    tool = "bamToFastq",
                    input = {"bam": target + "_100k.bam"},
                    output = {"fastq": target + "_100k.fastq.tmp"}))
        sample_reads(workflow, conf, 100000, "fastq")
        attach_back(workflow, PythonCommand(stat_fastqStat,
            input = {"seq": [ t + "_100k.seq" for t in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_seq_quality.json"},
            param = {"samples": conf.treatment_bases,
                     "seq_type": conf.seq_type.strip().split(",")[0].strip().lower()}))
        attach_back(workflow, PythonCommand(
            seq_quality_doc,
            input = {"tex": tex, "json": conf.json_prefix + "_seq_quality.json"},
            output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
            param = {"seq_type": conf.seq_type.strip().split(",")[0].strip().lower(), "reps": len(conf.treatment_pairs), "se_samples": conf.treatment_bases}))

    elif conf.seq_type.startswith("bed"): ## not evaluate sequence quality
        _versatile_format(workflow, conf)