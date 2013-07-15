__author__ = 'qinqianhappy'

#########################################################
#
# union DHS evaluation
#
#########################################################
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.helpers import *

def union_DHS_overlap(workflow, conf, tex):
    """ use hotspot d narrow peaks for DHS evaluation """
    for i, target in enumerate(conf.treatment_targets):
        if conf.peakcalltool == "hotspot":
            peaks = conf.hotspot_reps_final_5M_prefix[i] + ".fdr0.01.pks.bed"
            output = target + "_DHS_overlap_peaks_bed"
        elif conf.peakcalltool == "macs2":
            peaks = target + "_5M_macs2_velcro_non_overlap_peaks.bed.final"
            output = target + "_macs2_DHS_overlap_peaks_bed"

        attach_back(workflow,
            ShellCommand(
                "{tool} -wa -u  \
                -a {input[pks_spot_bed]} -b {input[DHS_peaks_bed]} > {output}",
                tool="intersectBed",
                input={"pks_spot_bed": peaks,
                       "DHS_peaks_bed": conf.get("lib", "dhs")},
                output = output,
                name = "Write out DHS overlap BED"))
    attach_back(workflow, PythonCommand(
        stat_dhs,
        input={"dhs_peaks": [ target + "_DHS_overlap_peaks_bed" if conf.peakcalltool == "hotspot" else target + "_macs2_DHS_overlap_peaks_bed"
                              for target in conf.treatment_targets ],
               "pks_spot_bed": [ conf.hotspot_reps_final_5M_prefix[i] + ".fdr0.01.pks.bed" if conf.peakcalltool == "hotspot" else target + "_5M_macs2_velcro_non_overlap_peaks.bed.final"
                                 for i, target in enumerate(conf.treatment_targets) ]},
        output = {"json": conf.json_prefix + "_dhs.json"},
        param= {"samples":conf.treatment_bases},
        name="DHS summary"))
    attach_back(workflow, PythonCommand(
        DHS_doc,
        input = {"tex": tex, "json": conf.json_prefix + "_dhs.json"},
        output = {"latex": conf.latex_prefix + "_dhs.tex"},
        param = {"reps": len(conf.treatment_pairs), "samples": conf.treatment_bases}))

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