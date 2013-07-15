#########################################################
#
# 1. library contamination QC
# 2. Reads mapping QC
#
#########################################################

import random, os
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from pkg_resources import resource_filename
from gcap.funcs.helpers import *

## bowtie command line for pair end and single end DNase
def _bowtie(workflow, conf):
    """
    Only all reads mapping, optionally
    """
    if conf.seq_type == "se":
        for raw, target in conf.treatment_pairs:
            bowtie = attach_back(workflow,
                ShellCommand(
                    "{tool} -n 1 -p {param[threads]} -S -m {param[max_align]} \
                    {param[genome_index]} {input[fastq]} {output[sam]} 2> {output[bowtie_summary]}",
                    input={"fastq": raw},
                    output={"sam": target + "_all.sam",
                            "bowtie_summary": target + "all_um_summary.txt"},
                    tool="bowtie",
                    param={"threads": 4,
                           "max_align": 1,
                           "genome_index": conf.get_path("lib", "genome_index")}))
            bowtie.update(param = conf.items("bowtie"))

    elif conf.seq_type == "pe":
        for n, target in enumerate(conf.treatment_targets):
            bowtie = attach_back(workflow,
                ShellCommand(
                    "{tool} -X 600 --chunkmbs 300 -n 1 -p {param[threads]} -S -m {param[max_align]} \
                    {param[genome_index]} -1 {input[fastq][0]} -2 {input[fastq][1]} {output[sam]} 2> {output[bowtie_summary]}",
                    input={"fastq": conf.treatment_raws[n]},
                    output={"sam": target + "_all.sam",
                            "bowtie_summary": target + "all_um_summary.txt"},
                    tool="bowtie",
                    param={"threads": 4,
                           "max_align": 1,
                           "genome_index": conf.get_path("lib", "genome_index")}))
            bowtie.update(param = conf.items("bowtie"))
    elif conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam") or conf.seq_type.startswith("bed"): ## skip mapping
        pass

## default bwa command line for pair end and single end DNase
def _bwa(workflow, conf):
    """bwa genome fasta genome index
      bwa aln genome.index.fa s_4_1_sequence.fastq > s_4_1.sai
      bwa aln genome.index.fa s_4_2_sequence.fastq > s_4_2.sai
      bwa sampe genome.index.fa s_4_1.sai s_4_2.sai s_4_1_sequence.fastq s_4_2_sequence.fastq > s_4.sam
    """
    if conf.seq_type == "pe":
        for n, target in enumerate(conf.treatment_targets):
            for pair in conf.treatment_raws[n]:
                attach_back(workflow, ShellCommand(
                    "{tool} aln -t {param[threads]} {input[index]} {input[fastq]} > {output[sai]}",
                    tool = "bwa",
                    input = {"index": conf.get("lib", "genome_index"),
                             "fastq": pair},
                    output = {"sai": pair + ".sai"},
                    param = {"threads": 4}))
            attach_back(workflow, ShellCommand(
                "{tool} sampe {input[index]} {input[sai][0]} {input[sai][1]} {input[fastq][0]} {input[fastq][1]} > {output[sam]}",
                tool = "bwa",
                input = {"index": conf.get("lib", "genome_index"),
                         "fastq": [ pair for pair in conf.treatment_raws[n] ],
                         "sai": [ pair + ".sai" for pair in conf.treatment_raws[n] ]},
                output = {"sam": target + "_all.sam"}))

    elif conf.seq_type == "se":
        for raw, target in conf.treatment_pairs:
            attach_back(workflow, ShellCommand(
                "{tool} aln -t {param[threads]} {input[index]} {input[fastq]} > {output[sai]}",
                tool = "bwa",
                input = {"index": conf.get("lib", "genome_index"),
                         "fastq": raw},
                output = {"sai": target + ".sai"},
                param = {"threads": 4}
            ))
            attach_back(workflow, ShellCommand(
                "{tool} samse {input[index]} {input[sai]} {input[fastq]} > {output[sam]}",
                tool = "bwa",
                input = {"index": conf.get("lib", "genome_index"),
                         "fastq": raw, "sai": target + ".sai"},
                output = {"sam": target + "_all.sam"}))

## parse bowtie summary for library contamination
def _bowtie_summary_parse(input=""):
    summary_content = open(input).read()
    print(summary_content)
    print("_" * 100)
    total_reads = int(re.findall(r"# reads processed: (.*)\n", summary_content)[0])

    # WARN: this is `unique mappable reads` as we set `-m 1`
    mappable_reads = int(re.findall(r"# reads with at least one reported alignment: (.*) \(", summary_content)[0])

    mappable_rate = float(mappable_reads) / float(total_reads)

    return {"total_reads": total_reads, "mappable_reads": mappable_reads, "mappable_rate": mappable_rate}

## summary of library contamination, if contaminated, exit
def stat_contamination(input = {"bowtie_summaries": [[]]},
                       output = {"json": ""}, param = {"samples": "", "species": "", "id": "", "correct_species": ""}):
    ## if library is contaminated, break
    library_contamination = {}
    library_contamination["header"] = {"sample": "sample name", "species": param["species"]}
    library_contamination["value"] = {}
    for a_summary, s in zip(input["bowtie_summaries"], map(underline_to_space, param["samples"])):
        ## each bowtie_summary has several species information
        library_contamination["value"][s] = {}
        for i, j in zip(a_summary, param["species"]):
            ## species 1, species2, species3
            parsed_summary = _bowtie_summary_parse(i)
            library_contamination["value"][s][j] = parsed_summary["mappable_rate"]

    for s in library_contamination["value"]:
        # s for different samples
        for i in param["species"]:
            # i for species, check each species
            if i != param["correct_species"]:
                if library_contamination["value"][s][param["correct_species"]] > library_contamination["value"][s][i] :
                    print("Library is correct !")
                else:
                    print("Library may be contaminated !")
                    sys.exit(1)

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict["stat"] = library_contamination
    json_dump(json_dict)

## library contamination QC, mapping to human, mouse and rat using bowtie for fast evaluation
def lib_contamination(workflow, conf, tex):
    """ input PE and SE sampled data 100 K to run,
     output json for library contamination, if not correct, shutdown
    """
    ## bowtie mapping back to mouse, rat and human
    if conf.seq_type == "se":
        for target in conf.treatment_targets:
            for species in dict(conf.items("contaminate")):
                attach_back(workflow,
                    ShellCommand(
                        "{tool} -p {param[threads]} -S -m {param[max_align]} \
                        {param[genome_index]} {input[fastq_sample]} {output[sam]} 2> {output[bowtie_summary]}",
                        input={"fastq_sample": target + "_100k.fastq"},
                        output={"sam": target + "_100k.sam",
                                "bowtie_summary": target + species + "_contam_summary.txt"},
                        tool="bowtie",
                        param={"threads": 4,
                               "max_align": 1,
                               "genome_index": conf.get_path("contaminate", species)}))
        all_species = [i for i, _ in conf.items("contaminate")]
        bowtie_summaries = []
        for target in conf.treatment_targets:
            bowtie_summaries.append([target + species + "_contam_summary.txt" for species in all_species])

        ## if library is contaminated, breaks
        attach_back(workflow,
            PythonCommand(stat_contamination,
                input={"bowtie_summaries": bowtie_summaries},
                output={"json": conf.json_prefix + "_contam.json"},
                param={"samples": conf.treatment_bases,
                       "id": conf.id,
                       "species": all_species,
                       "correct_species": conf.get("Basis", "species")}))

    elif conf.seq_type == "pe":
        for n, target in enumerate(conf.treatment_pair_data):
            for species in dict(conf.items("contaminate")):
                attach_back(workflow,
                    ShellCommand(
                        "{tool} -X 600 --chunkmbs 300 -n 1 -p {param[threads]} -S -m {param[max_align]} \
                        {param[genome_index]} -1 {input[fastq_sample][0]} -2 {input[fastq_sample][1]} {output[sam]} 2> {output[bowtie_summary]}",
                        input={"fastq_sample": [ t + "_100k.fastq" for t in target ]},
                        output={"sam": conf.treatment_targets[n] + "_100k.sam",
                                "bowtie_summary": conf.treatment_targets[n] + species + "_contam_summary.txt"},
                        tool="bowtie",
                        param={"threads": 4,
                               "max_align": 1,
                               "genome_index": conf.get_path("contaminate", species)},
                        name = "library contamination"))
        all_species = [i for i, _ in conf.items("contaminate")]
        bowtie_summaries = []
        for target in conf.treatment_targets:
            bowtie_summaries.append([target + species + "_contam_summary.txt" for species in all_species])

        ## if library is contaminated, breaks
        attach_back(workflow,
            PythonCommand(stat_contamination,
                input={"bowtie_summaries": bowtie_summaries},
                output={"json": conf.json_prefix + "_contam.json"},
                param={"samples": conf.treatment_bases,
                       "id": conf.id,
                       "species": all_species,
                       "correct_species": conf.get("Basis", "species")}))

## parse mapping json file to latex report
def reads_doc(input = {"tex": "", "json_autosome": ""},
              output = {"raw": "", "mapping": ""}, param = {"seq_type": "", "reps": "", "samples": "", "tool": ""}):
    """
    these values will be replaced by Jim's codes,
    currently by bowtie and fastqc,
    use autosome reads mappable ratio
    """
    json_au = json_load(input["json_autosome"])["stat"]  ## autosome
    ## json = json_load(input["json"])["stat"]
    ## use correct order

    if param["seq_type"] == "se":
        raw =  [ count_in_million(int(json_au[s]["total"])) for s in param["samples"] ]
        total = sum([int(json_au[j]["total"]) for j in json_au])
        mapped_number = [ int(json_au[s]["map"]) for s in param["samples"] ]
    else:
        raw =  [ "%s, %s" % (count_in_million(int(int(json_au[s]["total"])/2)),
                             count_in_million(int(int(json_au[s]["total"])/2))) for s in param["samples"] ]
        total = sum([int(json_au[j]["total"]) for j in json_au])
        mapped_number = [int(json_au[s]["map"]) for s in param["samples"]]
    ## use autosome ratio
    mapped_rate = [ decimal_to_latex_percent(json_au[s]["ratio"]) for s in param["samples"] ]
    ## raw reads doc
    reads_latex = JinjaTemplateCommand(
        template = input['tex'],
        name = "raw_reads",
        param = {"section_name": "raw_reads",
                 "render_dump": output["raw"],
                 "reads": raw,
                 "reps": param["reps"],
                 "combo": count_in_million(total)})
    template_dump(reads_latex)

    ## mapping status doc
    mapping_latex = JinjaTemplateCommand(
        template = input['tex'],
        name = "reads mapping",
        param = {"section_name": "mapping",
                 "render_dump": output["mapping"],
                 "map": mapped_rate, "combo": count_in_million(sum(mapped_number)),
                 "tool": param["tool"],
                 "reps": param["reps"]})
    template_dump(mapping_latex)

## write mapping information into json file
def autosome_map(input = {"count": ""}, output = {"json": ""}, param = {"samples": ""}):
    """
    use autosome mapping ratio to replace all reads mapping ratio
    """
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    for i, s in zip(input["count"], param["samples"]):
        json_dict["stat"][s] = {}
        with open(i, 'rU') as f:
            data = f.read().strip().split("\t")
            no_chrM_tag = data[0]
            total = data[1]

        json_dict["stat"][s]["ratio"] = round(float(no_chrM_tag) / float(total), 3)
        json_dict["stat"][s]["map"] = no_chrM_tag
        json_dict["stat"][s]["total"] = total

    json_dump(json_dict)

## reads mapping by customized mapping tool, calculate autosome(exclude M, X, Y) mapping ratio as we discussed before
def reads_mapping(workflow, conf, tex):
    """
    reads mapping for all reads
    and convert sam to bam, sort
    Sam or Bam input, only do statistics
    """
    if conf.maptool == "bowtie":
        _bowtie(workflow, conf)
    elif conf.maptool == "bwa":
        _bwa(workflow, conf)

    ## filter mitochondria, chrX, chrY tags
    if not conf.seq_type.startswith("bed"): ## reads BED files do not need to map
        if conf.maptool == "bowtie":
            uniq_tag = resource_filename("gcap", "pipeline-scripts/autosome_uniq_map_bowtie.awk")
        else:
            uniq_tag = resource_filename("gcap", "pipeline-scripts/autosome_uniq_map_bowtie.awk")
        for target in conf.treatment_targets:
            attach_back(workflow, ShellCommand(
                "{tool} -f {param[script]} {input[sam]} > {output[count]}",
                tool = "awk",
                input = {"sam": target + "_all.sam"},
                output = {"count": target + "_uniq_tag_autosome.count"},
                param = {"script": uniq_tag}))
        attach_back(workflow, PythonCommand(
            autosome_map,
            input = {"count": [ target + "_uniq_tag_autosome.count" for target in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_um_autosome.json"}, param = {"samples": conf.treatment_bases}))

        ## raw reads and reads mapping qc
        attach_back(workflow, PythonCommand(
            reads_doc,
            input = {"tex": tex, ## "json": conf.json_prefix + "_mapping.json", ## for uniform reads statistics
                     "json_autosome": conf.json_prefix + "_um_autosome.json"},
            output = {"raw": conf.latex_prefix + "raw.tex",
                      "mapping": conf.latex_prefix + "_mapping.tex"},
            param = {"seq_type": conf.seq_type,
                     "reps": len(conf.treatment_pairs),
                     "samples": conf.treatment_bases,
                     "tool": conf.maptool}))

        ## convert to BAM file
        for n, target in enumerate(conf.treatment_targets):
            if not conf.seq_type.startswith("bam"):
                attach_back(workflow,
                    ShellCommand(
                        "{tool} view -bt {input[chrom_len]} {input[sam]} -o {output[bam]}",
                        tool="samtools",
                        input={"sam": target + "_all.sam", "chrom_len": conf.get_path("lib", "chrom_len")},
                        output={"bam": target + ".bam"}))
                workflow.update(param=conf.items("lib"))
            else:
                ## BAM input Files do not need to convert
                link = attach_back(workflow, ShellCommand(
                    "{tool} -sf {input[bam]} {output[bam]}",
                    tool = "ln",
                    input = {"bam": conf.treatment_bam[n]},
                    output = {"bam": target + ".bam"},   ## actually, may not be sorted bam files
                    param = None))