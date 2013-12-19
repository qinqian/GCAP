__author__ = 'qinqianhappy'
#########################################################
#
# library complexity evaluation
#   1. SAM file to sorted BAM files
#   2. 5M library complexity by census
#
#########################################################
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.helpers import *
from gcap.funcs.sampling import *

def library_complexity(workflow, conf, tex):
    """  sampling 5M reads for following analysis
    use 5M raw reads for estimation
    """
    ## choose sampling modes
    if not conf.seq_type.startswith("bed"):
        ## top 5M raw reads, output SAM
        if conf.seq_type.startswith("bam"):
            sample_reads(workflow, conf, 5000000, "bam")
        else:
            sample_reads(workflow, conf, 5000000, "sam")
        ## for fastq files or sam input
        for target in conf.treatment_targets:
            lib_sort = attach_back(workflow, ShellCommand(
                "{tool} -Xmx5g -XX:ParallelGCThreads={param[threads]} -jar {param[sort]} I={input[bam]} O={output[bam]} SO=coordinate \
                VALIDATION_STRINGENCY=SILENT",
                tool = "java",
                input = {"bam": target + "_5M.sam"},
                output = {"bam": target + "_5M_sort.bam"},
                param = {"threads": 4},
                name = "Library Sort Sam to Bam"))
            lib_sort.update(param = conf.items("picard"))
        for target in conf.treatment_targets:
        #            lib_dup = attach_back(workflow, ShellCommand(
        #                "{tool} -Xmx5g -XX:ParallelGCThreads={param[threads]} -jar {param[markdup]} I={input[bam]} O={output[bam]} METRICS_FILE={output[metrics]} REMOVE_DUPLICATES=false \
        #                VALIDATION_STRINGENCY=SILENT",
        #                tool = "java",
        #                input = {"bam": target + "_picard_sort.bam.5M"},
        #                output = {"bam": target + "_markdup.bam.5M", "metrics": target + "_markdup_metric.5M"},
        #                param = {"threads": 4}))
        #            lib_dup.update(param = conf.items("picard"))
        #        attach_back(workflow, PythonCommand(
        #            stat_redun_picard,
        #            input = {"picard": [target + "_markdup_metric.5M" for target in conf.treatment_targets ]},
        #            output = {"json": conf.json_prefix + "_redun.json"}, param = {"samples": conf.treatment_bases, "format": "notbed"}))
            ## use census to estimation instead of picard
            attach_back(workflow, ShellCommand(
                "{tool} {param[hist]} {param[se_or_pe]} {param[exclude]} {input[sorted_bam]} | {tool} {param[calc]} - > {output[metrics]}",
                tool = "python2.7",
                input = {"sorted_bam": target + "_5M_sort.bam"},
                output = {"metrics": target + "_census.metric"},
                param = {"exclude": conf.get("census", "census_exclude"),  ## hg19 exclude only now, add mouse exclude latter
                         "hist": conf.get("census", "hist"),
                         "calc": conf.get("census", "calc"),
                         "se_or_pe": "-s" if "se" in conf.seq_type else " "},
                name = "census program"))

#            ## PBC1, PBC2
#            attach_back(workflow, ShellCommand(
#                '''file={input[histo]}
#                work=$(basename $file)
#                prefix=${work%.*}
#                suffix=$3
#                proj=${2}/${prefix}
#                {tool} '{l+=1;ul[$1"\t"$2"\t"$3"\t"$6]+=1} END {for (u in ul) print ul[u]}' ${proj}.${suffix} |\
#                {tool} '{n[$1]+=1} END {for (i in n) print i"\t"n[i]}' > ${proj}.${suffix}.histo
#                ''',
#                tool = "awk", name = "PBC"))
#
#            ## Preseq
#            attach_back(workflow, ShellCommand(
#                '''
#                sort -k1n ${proj}.${suffix}.histo > ${proj}.${suffix}.histo.sort
#                ## time c_curve -o ${proj}.lib_complexity -H ${proj}.${suffix}.histo.sort
#                ## library yield, extract up, so slow, not calculated yet
#                ## lc_extrap -b 90 -e 10000000 -s 1000000 ${proj}.sort.bed -o ${proj}.future_yield
#                time {tool} -H ${proj}.${suffix}.histo.sort -o ${proj}.future_yield
#                ''',
#                tool = "lc_extrap"))

        attach_back(workflow, PythonCommand(
            stat_redun_census,
            input = {"census": [target + "_census.metric" for target in conf.treatment_targets]},
            output = {"json": conf.json_prefix + "_redun.json"},
            param = {"samples": conf.treatment_bases, "format": "notbed"},
            name = "census stat"))

        attach_back(workflow, PythonCommand(
            redundancy_doc,
            input = {"tex": tex, "json": conf.json_prefix + "_redun.json"},
            output = {"redun": conf.latex_prefix + "redun.tex"},
            param = {"reps": len(conf.treatment_pairs), "samples": conf.treatment_bases},
            name = "census document"))

    else:  ## get top 5M reads, use shell to calculate redundancy simply
        sample_reads(workflow, conf, 5000000, "bed")

## parse picard results to json
def stat_redun_picard(input = {"picard": "", "built_count": ""}, output = {"json": ""}, param = {"samples": "", "format": ""}):
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
            total_loc = open(b[0]).read().strip().split()[0]
            uniq_loc = open(b[1]).read().strip().split()[0]
            redun_ratio = 1 - float(uniq_loc) / float(total_loc)
            json_dict["stat"][s] = redun_ratio
    json_dump(json_dict)

## parse census results to json
def stat_redun_census(input = {"census": "", "redundancy": ""}, output = {"json": ""}, param = {"samples": "", "format": ""}):
    """
    '\tTotal unique reads (molecules):\t2675277 (99.8%)\n'
    1 - 99.8%
    this is for BED awk statistics, too
    """
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    if param["format"] != "bed":
        for f, s in zip(input["census"], param["samples"]):
            data = open(f).readlines()[4]
            d = 1 - float(re.findall("\((\d+.\d+)\%\)", data, re.I)[0])/100
            json_dict["stat"][s] = d
    else:
        ## for BED reads input, use awk to count
        for b, s in zip(input["redundancy"], param["samples"]):
            redun_ratio = open(b).read().strip().split()[0]
            json_dict["stat"][s] = redun_ratio
    json_dump(json_dict)

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
