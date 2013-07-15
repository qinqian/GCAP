
#########################################################
#
# 1. sort peaks by score
# 2. peaks conservation in 3 gradients
#
#########################################################
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.helpers import *

def conservation(workflow, conf, tex):
    ## use 5M reads hotspot b for evaluation of conservation
    for target in conf.treatment_targets:
        ## get non-promotor peaks regions
        if conf.peakcalltool == "hotspot":
            peaks = target + "_5M_velcro_non_overlap_hotspot.bed.final"
            non_promotor = target + "_5M_non-promotor_peaks_bed"
            phascon = target + "_100bp_phastcon.score"
        elif conf.peakcalltool == "macs2":
            peaks = target + "_5M_macs2_velcro_non_overlap_peaks.bed.final"
            non_promotor = target + "_macs2_5M_non-promotor_peaks_bed"
            phascon = target + "_macs2_100bp_phastcon.score"
        attach_back(workflow,
            ShellCommand(
                "{tool} -v -a {input[peaks]} -b {input[tss]} > {output}",
                tool = "intersectBed",
                input = {"peaks": peaks, "tss": conf.get("lib", "tss")},
                output = non_promotor,
                name = "Get non-promotor peaks regions",
                param = None))

        get_top_peaks = attach_back(workflow,
            ShellCommand(
                "{tool} -r -g -k 5 {input} | head -n {param[peaks]} > {output}",
                tool="sort",
                input=non_promotor,
                output=non_promotor + ".top1000",
                param={'peaks': 1000}, name="top summits for conservation"))
        get_top_peaks.update(param=conf.items('conservation'))

        ## get summits is implemented in conservation_average.py
        attach_back(workflow,
            ShellCommand(
                "{tool} -w {param[width]} -d {input[phastcon_db]} {input[bed]} 1>{output}",
                tool = "conservation_average.py",
                input= {"phastcon_db": conf.get_path("lib", "phast"),
                        "bed": non_promotor + ".top1000"},
                output = phascon,
                param = {"width": 100}))

    ## Phastcon json
    attach_back(workflow,
        PythonCommand(stat_conserv,
            input = {"phastcon":[ target + "_100bp_phastcon.score" if conf.peakcalltool == "hotspot" else target + "_macs2_100bp_phastcon.score"
                                  for target in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_conserv.json"},
            param = {"sample": conf.treatment_bases}))

    ## Phaston latex
    attach_back(workflow,
        PythonCommand(conserv_doc,
            input = {"tex": tex, "json": conf.json_prefix + "_conserv.json"},
            output = {"latex": conf.latex_prefix + "_conserv.tex"},
            param = {"reps": len(conf.treatment_pairs),
                     "sample": conf.treatment_bases}))


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