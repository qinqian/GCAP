__author__ = 'qinqianhappy'

##################################################################
#
# replicates consistency
#  1. I/U intersect BED over union BED(filtered hotspot regions) by bedtools
#  2. bigwiggle correlation on filtered DHS regions or genome wide
#
###################################################################

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.helpers import *

def peaks_reps_preprocess(workflow, conf):
    """
    remove blacklist and outlier for hotspot b
    """
    has_velcro = conf.get("lib", "velcro")  ## mouse genome does not have blacklist available
    for i, target in enumerate(conf.treatment_targets):
        ## remove blacklist for human
        if conf.peakcalltool == "hotspot":
            input_peaks = conf.hotspot_reps_final_5M_prefix[i] + ".hot.bed"
            non_velcro = target + "_5M_velcro_non_overlap_hotspot.bed"
        elif conf.peakcalltool == "macs2":
            input_peaks = target + "_5M_macs2_peaks.bed"
            non_velcro = target + "_5M_macs2_velcro_non_overlap_peaks.bed"

        ## awk and bedClip to remove outlier location from above input
        attach_back(workflow,
            ShellCommand(
                "sed 1d {input} | {tool} '{{if ($2 >= 0 && $2 < $3) print}}' - | sort-bed -  | bedops -m - > {output}",
                tool="awk",
                input = input_peaks,
                output = input_peaks + ".tmp",
                name = "filter bed files outlier location"))

        # prototype used here to do the similar thing on bedclip
        bed_clip = attach_back(workflow,
            ShellCommand(
                template="{tool} {input} {param[chrom_len]} {output}",
                tool="bedClip",
                input=input_peaks + ".tmp",
                output = input_peaks + ".final",
                param={'chrom_len': conf.get_path("lib", "chrom_len")},
                name="bedclip filter"))
        bed_clip.allow_fail = True

        if has_velcro:
            attach_back(workflow,
                ShellCommand(
                    "{tool} -v -a {input[peaks]} -b {input[velcro_peaks_bed]} > {output}",
                    tool="intersectBed",
                    input={"peaks": input_peaks + ".final",
                           "velcro_peaks_bed": conf.get("lib", "velcro")},
                    output = non_velcro + ".final",
                    name = "filter blacklist",
                    param=None))
        else: ## for mouse, just rename to *_5M_velcro_non_overlap_peaks.bed
            attach_back(workflow,
                ShellCommand(
                    "cp {input[peaks]} {output}",
                    tool="cp",
                    input={"peaks": input_peaks + ".final"},
                    output = non_velcro + ".final",
                    name = "No blacklist, not filter blacklist"))

def stat_reps(input={"5M_overlap": "", "5M_cor": "", "union": ""},
              output={"json": ""}, param={"cor": "genome"}):
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
            if param["cor"] == "genome":
                rep_cor = f.read().strip().split()[2]
            else:
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

def peaks_reps_evaluating(workflow, conf, tex):
    ## use 5M reads hotspot filtered regions for evaluation, intersect region ratio in pairwise ways to calculate I/U
    for i in range(len(conf.treatment_targets)):
        for j in range(i+1, len(conf.treatment_targets)):
            if conf.peakcalltool == "hotspot":
                one_rep = conf.treatment_targets[i] + "_5M_velcro_non_overlap_hotspot.bed.final"
                another_rep = conf.treatment_targets[j] + "_5M_velcro_non_overlap_hotspot.bed.final"
            elif conf.peakcalltool == "macs2":
                one_rep = conf.treatment_targets[i] + "_5M_macs2_velcro_non_overlap_peaks.bed.final"
                another_rep =  conf.treatment_targets[j] + "_5M_macs2_velcro_non_overlap_peaks.bed.final"
            attach_back(workflow,
                ShellCommand(
                    "{tool} -a {input[one_rep]} -b {input[another_rep]} > {output}",
                    tool="intersectBed",
                    input = {"one_rep": one_rep,
                             "another_rep": another_rep},
                    output= conf.treatment_targets[i] + "_reps_overlap_" + "rep" + str(i+1) + "rep" + str(j+1) + ".bed",
                    name = "replicates Overlap"))
    ## Implemented by Jim
    ## correlation in prepared union DHS regions(filtered by blacklist), bigWigCorrelate -restrict=testid_reps_union_region.bb testid_treat_rep1_5M.bw testid_treat_rep2_5M.bw
    ## correlation in genome wide, wigCorrelate,  wigCorrelate one.wig two.wig ... n.wig
    for i in range(len(conf.treatment_targets)):
        for j in range(i+1, len(conf.treatment_targets)):
            if conf.peakcalltool == "macs2":
                one_rep = conf.treatment_targets[i] + "_5M_macs2_treat.bw"
                another_rep = conf.treatment_targets[j] + "_5M_macs2_treat.bw"
            elif conf.peakcalltool == "hotspot":
                one_rep = conf.treatment_targets[i] + "_5M.bw"
                another_rep = conf.treatment_targets[j] + "_5M.bw"

            if conf.get("tool", "cor").strip().lower() == "genome":
                attach_back(workflow,
                    ShellCommand(
                        "{tool} {input[one_rep]} {input[another_rep]} 1>{output[cor_score]}",
                        tool="wigCorrelate",
                        input = {"one_rep": one_rep,
                                 "another_rep": another_rep,
                                 "union": conf.get("lib", "filtered_dhs_bb")},
                        output= {"cor_score": conf.treatment_targets[i] + "_reps_" + str(i+1) + str(j+1) + "cor_score"},
                        name = "replicates bigwiggle genome wide correlation"))

            elif conf.get("tool", "cor").strip().lower() == "union":
                attach_back(workflow,
                    ShellCommand(
                        "{tool} -restrict={input[union]} {input[one_rep]} {input[another_rep]} 1>{output[cor_score]}",
                        tool="bigWigCorrelate",
                        input = {"one_rep": one_rep,
                                 "another_rep": another_rep,
                                 "union": conf.get("lib", "filtered_dhs_bb")},
                        output= {"cor_score": conf.treatment_targets[i] + "_reps_" + str(i+1) + str(j+1) + "cor_score"},
                        name = "replicates bigwiggle union DHS correlation"))
                ## peaks replicates json
    bed = [ (conf.treatment_targets[i] + "_reps_overlap_" + "rep" + str(i+1) + "rep" + str(j+1) + ".bed", "Rep %s %s" %(i + 1, j + 1))
            for i in range(len(conf.treatment_targets)) for j in range(i+1, len(conf.treatment_targets)) ]

    cor = [ (conf.treatment_targets[i] + "_reps_" + str(i+1) + str(j+1) + "cor_score", "Rep %s %s" %(i + 1, j + 1))
            for i in range(len(conf.treatment_targets)) for j in range(i+1, len(conf.treatment_targets)) ]

    ## get 5M union hotspot regions for regions number I/U jaccard index evaluation
    union = attach_back(workflow, ShellCommand(
        "{tool} -m {param[reps]} > {output[union]}",
        tool = "bedops",
        input = {"beds": [ target + "_5M_velcro_non_overlap_hotspot.bed.final" if conf.peakcalltool == "hotspot" else target + "_5M_macs2_velcro_non_overlap_peaks.bed.final"
                           for target in conf.treatment_targets ]},
        output = {"union": conf.prefix + "_hotspot_5M_union.bed"}))
    union.param.update({"reps": " ".join(union.input['beds'])})

    ## get output to json
    attach_back(workflow, PythonCommand(
        stat_reps,
        input = {"5M_overlap": bed, "5M_cor": cor, "union": conf.prefix + "_hotspot_5M_union.bed"},
        output = {"json": conf.json_prefix + "_reps.json"},
        param = {"cor": conf.get("tool", "cor")}))

    ## latex document
    attach_back(workflow, PythonCommand(
        reps_doc,
        input = {"json": conf.json_prefix + "_reps.json", "tex": tex},
        output = {"latex": conf.latex_prefix + "_reps.tex"},
        param = {"samples": conf.treatment_bases}))
