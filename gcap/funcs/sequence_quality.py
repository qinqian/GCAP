#############################################################
# The first step for GCAP QC, add multiple format support
#    input: fastq, use fastqStatsAndSubsample
#    output: sequence quality and subsampled reads
#############################################################

from samflow.command import ShellCommand, PythonCommand
from gcap.funcs.helpers import *


def read_quality(workflow, conf, tex):
    if conf.pe:
        for raw, target in conf.treatment_pairs_pe:
            attach_back(workflow,
                        ShellCommand(
                            "{tool} {input[fastq][0]} {input[fastq][1]} {output[stat][0]} {output[stat][1]",
                            tool = "dac_pe_read_quality",
                            input = {"fastq": raw},
                            output = {"stat": [ i + "_read_quality.qc" for i in target ]}))

        # attach_back(workflow, PythonCommand(stat_fastqStat,
        #                                     input = {"seq": [ [ p + "_100k.seq" for p in target ] for target in conf.treatment_pair_data ]},
        #                                     output = {"json": conf.json_prefix + "_seq_quality.json"},
        #                                     param = {"samples": conf.treatment_bases, "seq_type": conf.pe}))
        # attach_back(workflow, PythonCommand(
        #     seq_quality_doc,
        #     input = {"tex": tex, "json": conf.json_prefix + "_seq_quality.json"},
        #     output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
        #     param = {"seq_type": conf.seq_type, "reps": len(conf.treatment_pairs),
        #              "pe_samples": conf.treatment_bases}))
    else:
        for raw, target in conf.treatment_pairs:
            sample_fq = {"stat": target + "_read_quality.qc"}
            attach_back(workflow,
                        ShellCommand(
                            "{tool} {input} {output[stat]}",
                            tool = "dac_se_read_quality",
                            input = raw,
                            output = sample_fq,
                            name = "100k read sequence quality and sequence length"))

        # attach_back(workflow, PythonCommand(stat_fastqStat,
        #                                     input = {"seq": [ t + "_100k.seq" for t in conf.treatment_targets ]},
        #                                     output = {"json": conf.json_prefix + "_seq_quality.json"},
        #                                     param = {"samples": conf.treatment_bases,
        #                                              "seq_type": conf.pe}))
        # attach_back(workflow, PythonCommand(
        #     seq_quality_doc,
        #     input = {"tex": tex, "json": conf.json_prefix + "_seq_quality.json"},
        #     output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
        #     param = {"seq_type": conf.seq_type, "reps": len(conf.treatment_pairs), "se_samples": conf.treatment_bases}))


## parse fastqStatsAndSubsample result
def stat_fastqStat(input = {"seq": ""}, output = {"json": ""}, param = {"samples": "", "seq_type": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    if param["seq_type"] == "pe":
        for i, s in zip(input['seq'], param['samples']):
            with open(i[0]) as f:
                data = f.readlines()
                json_dict["stat"][s + "_pair1"] = {}
                json_dict["stat"][s + "_pair1"]['count'] = data[0].strip().split()[1]
                json_dict["stat"][s + "_pair1"]['basecount'] = data[1].strip().split()[1]
                json_dict["stat"][s + "_pair1"]['quality'] = round(float(data[8].strip().split()[1]), 1)
                json_dict["stat"][s + "_pair1"]['std'] = data[9].strip().split()[1]
                json_dict["stat"][s + "_pair1"]['len'] = data[4].strip().split()[1]
            with open(i[1]) as f:
                data = f.readlines()
                json_dict["stat"][s + "_pair2"] = {}
                json_dict["stat"][s + "_pair2"]['count'] = data[0].strip().split()[1]
                json_dict["stat"][s + "_pair2"]['basecount'] = data[1].strip().split()[1]
                json_dict["stat"][s + "_pair2"]['quality'] = round(float(data[8].strip().split()[1]), 1)
                json_dict["stat"][s + "_pair2"]['std'] = data[9].strip().split()[1]
                json_dict["stat"][s + "_pair2"]['len'] = data[4].strip().split()[1]
    elif param["seq_type"] == "se":
        for i, s in zip(input['seq'], param['samples']):
            with open(i) as f:
                data = f.readlines()
                json_dict["stat"][s] = {}
                json_dict["stat"][s]['count'] = data[0].strip().split()[1]
                json_dict["stat"][s]['basecount'] = data[1].strip().split()[1]
                json_dict["stat"][s]['quality'] = data[8].strip().split()[1]
                json_dict["stat"][s]['std'] = data[9].strip().split()[1]
                json_dict["stat"][s]['len'] = data[4].strip().split()[1]
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
