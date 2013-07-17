#########################################################
#
# fragment size evaluation
#  1. pair end picard
#  2. single end macs2
#########################################################

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.sampling import *
from gcap.funcs.helpers import *

## fragment picard and macs2 command line
def fragment_size(workflow, conf, tex):
    """
    picard for PE,
    MACS2 for SE, because of background noise, we choose alternative fragment size closest to size selection
    """
    ## picard for pair end
    if conf.seq_type == "pe":
        for target in conf.treatment_targets:
            insert_size = attach_back(workflow, ShellCommand(
                "{tool} -Xmx5g -XX:ParallelGCThreads={param[threads]} -jar {param[insertsize]} HISTOGRAM_FILE={output[pdf]} I={input[bam]} O={output[insert]} VALIDATION_STRINGENCY=SILENT",
                tool = "java",
                input = {"bam": target + "_5M_sort.bam"},
                output = {"pdf": target + "_picard_insert_5M.pdf", "insert": target + "_insert_metric.5M"},
                param = {"threads": 4},
                name = "Fragment size picard"
            ))
            insert_size.update(param = conf.items("picard"))

        attach_back(workflow, PythonCommand(
            stat_frag_std,
            input = {"insert": [ target + "_insert_metric.5M" for target in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_frag.json"},
            param = {"samples": conf.treatment_bases,
                     "frag_tool": "picard"}))

    elif conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam"):
        if conf.seq_type.split(",")[1].strip() == "pe":
            for target in conf.treatment_targets:
                insert_size = attach_back(workflow, ShellCommand(
                    "{tool} -Xmx5g -XX:ParallelGCThreads={param[threads]} -jar {param[insertsize]} HISTOGRAM_FILE={output[pdf]} I={input[bam]} O={output[insert]} VALIDATION_STRINGENCY=SILENT",
                    tool = "java",
                    input = {"bam": target + "_5M_sort.bam"},
                    output = {"pdf": target + "_picard_insert_5M.pdf", "insert": target + "_insert_metric.5M"},
                    param = {"threads": 4}))
                insert_size.update(param = conf.items("picard"))
            attach_back(workflow, PythonCommand(
                stat_frag_std,
                input = {"insert": [ target + "_insert_metric.5M" for target in conf.treatment_targets ]},
                output = {"json": conf.json_prefix + "_frag.json"},
                param = {"samples": conf.treatment_bases,
                         "frag_tool": "picard"}))
        elif conf.seq_type.split(",")[1].strip() == "se":
            for target in conf.treatment_targets:
                fragment_size = attach_back(workflow, ShellCommand(
                    "{tool} predictd -i {input[bam]} --rfile {param[prefix]}",
                    tool = "macs2",
                    input = {"bam": target + "_5M_sort.bam"},
                    output = {"R": target + "_5M_model.R"},
                    param = {"prefix": target + "_5M"}))
            ## extract standard deviation from MACS2 model.R,
            ## use m, p, and pileup value for standard deviation; mean fragment size is provided(choose the one with highest correlation)
            attach_back(workflow, PythonCommand(
                stat_frag_std,
                input = {"r": [ target + "_5M_model.R" for target in conf.treatment_targets ]},
                output = {"json": conf.json_prefix + "_frag.json", "r": [ target + "_5M_frag_sd.R" for target in conf.treatment_targets ] },
                param = {"samples": conf.treatment_bases,
                         "frag_tool": "macs2"},
                name = "macs2 model R script parser"))
    elif conf.seq_type.startswith("bed"):
        for target in conf.treatment_targets:
            fragment_size = attach_back(workflow, ShellCommand(
                "{tool} predictd -i {input[bed]} --rfile {param[prefix]}",
                tool = "macs2",
                input = {"bed": target + "_5M.bed"},
                output = {"R": target + "_5M_model.R"},
                param = {"prefix": target + "_5M"}))
            attach_back(workflow, PythonCommand(
                stat_frag_std,
                input = {"r": [ target + "_5M_model.R" for target in conf.treatment_targets ]},
                output = {"json": conf.json_prefix + "_frag.json", "r": [ target + "_5M_frag_sd.R" for target in conf.treatment_targets ] },
                param = {"samples": conf.treatment_bases,
                         "frag_tool": "macs2"}))

    elif conf.seq_type == "se":  ## for single end fastq input and SE, PE bam input, use macs2
        for target in conf.treatment_targets:
            fragment_size = attach_back(workflow, ShellCommand(
                "{tool} predictd -i {input[bam]} --rfile {param[prefix]}",
                tool = "macs2",
                input = {"bam": target + "_5M_sort.bam"},
                output = {"R": target + "_5M_model.R"},
                param = {"prefix": target + "_5M"}))
        attach_back(workflow, PythonCommand(
            stat_frag_std,
            input = {"r": [ target + "_5M_model.R" for target in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_frag.json", "r": [ target + "_5M_frag_sd.R" for target in conf.treatment_targets ] },
            param = {"samples": conf.treatment_bases,
                     "frag_tool": "macs2"}))
    attach_back(workflow, PythonCommand(
        frag_doc,
        input = {"json": conf.json_prefix + "_frag.json", "tex": tex},
        output = {"latex": conf.latex_prefix + "_frag.tex"},
        param = {"reps": len(conf.treatment_pairs),
                 "samples": conf.treatment_bases}))

## parse macs2 model R script
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
            if line.startswith("xcorr"):
                values['xcorr'] = line
            if line.startswith("ycorr"):
                values['ycorr'] = line
    return values

## calculate fragment mean length and standard variation
def stat_frag_std(input = {"r": "", "insert": ""}, output = {"json": "", "r": ""}, param = {"samples": "", "frag_tool": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    if param["frag_tool"] == "macs2":
        ## single end
        for rin, rout, s in zip(input["r"], output["r"], param["samples"]):
            values = get_size(rin)
            with open(rout, 'w') as f:
                f.write(values['positive'])
                f.write(values['minus'])
                f.write(values['xcorr'])
                f.write(values['ycorr'])
                f.write("xcorr.max = xcorr[which(ycorr==max(ycorr))]\n")
                f.write(values['x'])
                f.write("p.expect = sum(x * p/100) \n")
                f.write("m.expect = sum(x * m/100) \n")
                f.write("p.sd = sqrt(sum(((x-p.expect)^2)*p/100)) \n")
                f.write("m.sd = sqrt(sum(((x-m.expect)^2)*m/100)) \n")
                f.write("cat(paste((p.sd + m.sd)/2, '\t', xcorr.max)) \n")

            std_frag = os.popen("Rscript %s" % rout).read().strip().split()
            json_dict["stat"][s] = "mean %s, sd %s" % (int(float(std_frag[1])), int(float(std_frag[0])))

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