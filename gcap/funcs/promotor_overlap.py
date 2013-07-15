__author__ = 'qinqianhappy'
#########################################################
#
# promotor overlap percentage
#
#########################################################

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.helpers import *

def _promoter_overlap(workflow, conf, tex):
    """ use bedtools to get promotor peaks enrichment,
    Genome Promotor: cutoff 50% promotor overlaps with mappability regions bases number
    Peaks Promotor: cutoff 50% promotor overlaps with peaks regions bases number
    use 5M reads hotspot b to evaluate promotor overlap
    """
    ## use hotspot mappable_region for both hotspot and macs2 peak calling
    attach_back(workflow,
        ShellCommand(
            "{tool} -e -{param[percentage]} {input[ref]} {input[map]} | wc -l > {output[count]}",
            tool = "bedops",
            input = {"ref":conf.get("lib", "tss"), "map": conf.get("hotspot", "mappable_region")},
            output = {"count": conf.prefix + "_promotor.count"},
            param = {"percentage": "50%"}))

    ## peaks region evaluation
    for i, target in enumerate(conf.treatment_targets):
        if conf.peakcalltool == "hotspot":
            peaks = target + "_5M_velcro_non_overlap_hotspot.bed.final"
        elif conf.peakcalltool == "macs2":
            peaks = target + "_5M_macs2_velcro_non_overlap_peaks.bed.final"
        attach_back(workflow,
            ShellCommand(
                "{tool} -e -{param[percentage]} {input[ref]} {input[promotor]} | wc -l > {output[count]} && wc -l {input[ref]} > {output[peaks]}",
                tool = "bedops",
                input = {"ref": peaks,
                         "promotor": conf.get("lib", "tss")},
                output = {"count": target + "_5M_promotor_overlap.count",
                          "peaks": target + "_peaks.count"},
                param = {"percentage": "50%"},
                name = "promotor overlaps"))
        ## promotor json
    attach_back(workflow, PythonCommand(
        stat_promotor,
        input = {"peaks_promotor": [ target + "_5M_promotor_overlap.count" for target in conf.treatment_targets ],
                 "peaks": [ target + "_peaks.count" for target in conf.treatment_targets ],
                 "promotor": conf.prefix + "_promotor.count",
                 "mappable": conf.get("hotspot", "mappable_region")},
        output = {"json": conf.json_prefix + "_promotor.json"},
        param = {"samples": conf.treatment_bases}))

    ## promotor latex document
    attach_back(workflow, PythonCommand(
        promotor_doc,
        input = {"json": conf.json_prefix + "_promotor.json",
                 "tex": tex},
        output = {"latex": conf.latex_prefix + "_promotor.tex"},
        param = {"samples": conf.treatment_bases,
                 "reps": len(conf.treatment_pairs)}))


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