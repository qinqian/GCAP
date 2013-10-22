__author__ = 'qinqianhappy'
#########################################################
#
# 1. NSC
# 2. RSC
#
#########################################################

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.sampling import *
from gcap.funcs.helpers import *

def strand_cor(workflow, conf, tex):
    """ use ccQualityControl run_spp.R to estimate
    Rscript run_spp.R -rf -c=<tagAlign/BAMfile> -savp -out=<outFile>
    support BAM/SAM, fastq input
    """
    for target in conf.treatment_targets:
        spp = attach_back(workflow, ShellCommand(
            "{tool} {param[spp]} -x={param[exclude]} -s={param[s]} -c={input[bam]} -out={output}",      ## no need for -savp, plot
            tool = "Rscript",
            input = {"bam": target + "_5M_sort.bam"},
            output = target + "_strand_cor",
            param = {"spp": conf.get("tool", "spp"),
                     "s": "-500:5:1500",
                     "exclude": "-500:-1"}))
        spp.param.update(conf.items("spp"))

    attach_back(workflow, PythonCommand(
        stat_strand_cor,
        input = {"metric": [target + "_strand_cor" for target in conf.treatment_targets]},
        output = {"json": conf.json_prefix + "_strand_cor.json"},
        param = {"samples": conf.treatment_bases}))

    attach_back(workflow, PythonCommand(
        strand_cor_doc,
        input = {"json": conf.json_prefix + "_strand_cor.json", "tex": tex},
        output = {"nsc_latex": conf.latex_prefix + "_nsc.tex", "rsc_latex": conf.latex_prefix + "_rsc.tex"},
        param = {"samples": conf.treatment_bases,
                 "reps": len(conf.treatment_pairs)}))

def stat_strand_cor(input = {"metric": ""}, output = {"json": ""}, param= {"samples": ""}):
    """RSC NSC: COL9 COL10, COL3: estFragLen"""
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}

    for i, j in zip(input["metric"], param["samples"]):
        data = open(i).readlines()
        l = data[0].strip().split("\t")
        json_dict["stat"][j] = {}
        json_dict["stat"][j]["RSC"] = l[8]
        json_dict["stat"][j]["NSC"] = l[9]
        json_dict["stat"][j]["Frag"] = l[2]
    json_dump(json_dict)

def strand_cor_doc(input = {"json": ""}, output = {"nsc_latex": "", "rsc_latex": ""}, param = {"samples": "", "reps": ""}):
    data = json_load(input["json"])["stat"]
    NSC = []
    RSC = []
    for i in param["samples"]:
        NSC.append(round(float(data[i]["NSC"]), 2))
        RSC.append(round(float(data[i]["RSC"]), 2))

    NSC_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "strand correlation latex",
        param = {"section_name": "NSC",
                 "render_dump": output["nsc_latex"],
                 "NSC": NSC,
                 "reps": param["reps"]})
    template_dump(NSC_latex)

    RSC_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "strand correlation latex",
        param = {"section_name": "RSC",
                 "render_dump": output["rsc_latex"],
                 "RSC": RSC,
                 "reps": param["reps"]})
    template_dump(RSC_latex)