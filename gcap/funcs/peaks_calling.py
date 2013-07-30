
__author__ = 'qinqianhappy'
###############################################################
#
# 1. peaks calling by hotspot v4 default, alternative macs2
# 2. generate 20bp resolution bigwiggle for all reads and 5M
# 3. generate peaks bed for all reads and 5M
# 4. SPOT score evaluation
# 5. peaks number evaluation
# 6. replicates consistency by calling reps_quality.py
# 7. calculate naive library complexity by awk for BED reads input
#
###############################################################
from gcap.funcs.reps_quality import peaks_reps_preprocess, peaks_reps_evaluating
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.helpers import *
from gcap.funcs.library_complexity import stat_redun_census, redundancy_doc

from pkg_resources import resource_filename

## for hotspot conf file and pipeline shell scripts
token_file = resource_filename("gcap", "static/runall.tokens.txt")
pipeline_scripts = resource_filename("gcap", "pipeline-scripts")

## call peaks by hotspot or macs2 to QC DNase data
def call_peaks(workflow, conf, tex):
    have_treat_reps = len(conf.treatment_pairs) >= 2 ## replicates
    if conf.peakcalltool == "hotspot":
        _hotspot_on_replicates(workflow, conf, tex)
    elif conf.peakcalltool == "macs2":
        _macs2_on_reps(workflow, conf, tex)
    ## filter replicates hotspot regions or macs2 peaks regions
    peaks_reps_preprocess(workflow, conf)
    if have_treat_reps:
        ## evaluate replicates consistency
        peaks_reps_evaluating(workflow, conf, tex)
        ## peaks calling on merged reads files
        if conf.peakcalltool == "hotspot":
            _hotspot_combo(workflow, conf)
        elif conf.peakcalltool == "macs2":
            _macs2_on_combo(workflow, conf)
    ## render peaks information to latex
    _peaks_calling_latex(workflow, conf, tex)

## make config files for hotspot program
def spot_conf(input = {"tag": "", "mappable_region": "", "spot_conf": "", "chrom_info": ""},
              output = {"dir": "", "conf": ""},
              param = {"K": "", "FDRS": "0.01", "species": "", "keep_dup": "T", "omit": ""}):
    cf = ConfigParser()
    ## preserve uppercase
    cf.optionxform = str
    cf.read(input["spot_conf"])
    setting = lambda cf, before, after: cf.set("script-tokenizer", before, after)
    setting(cf,"_TAGS_", input["tag"])
    setting(cf,"_GENOME_", param["species"])
    setting(cf, "_K_", param["K"])
    setting(cf,"_MAPPABLE_FILE_", input["mappable_region"])
    setting(cf, "_CHROM_FILE_", input["chrom_info"])
    setting(cf, "_DENS_", output["dir"])
    setting(cf, "_OUTDIR_", output["dir"])
    setting(cf, "_CLEAN_", "T")                       ## clean up output directory
    setting(cf, "_RANDIR_", output["dir"])
    setting(cf, "_FDRS_", param["fdrs"])
    setting(cf, "_DUPOK_", param["keep_dup"])         ## keep duplicates
    setting(cf, "_CHKCHR_", "chrX")     ## To check if data contains this chromosome
    setting(cf, "_CHECK_", "F")         ## change to F in case that hotspot error for no wavelets peaks, other method
    setting(cf, "_OMIT_REGIONS_", param["omit"])
    cf.write(open(output["conf"], "w"))

## collect spot 5M reads for each replicates into json files
def stat_spot_on_replicates(input = {"spot_files": ""}, output = {"json": ""},
                            param = {"samples": ""}):
    spot = {}
    spot["input"] = input
    spot["output"] = output
    spot["stat"] = {}
    for spot_file, s in zip(input["spot_files"], param["samples"]):
        with open(spot_file, "rU") as f:
            stat = f.readlines()[1].strip().split()
            stat_dict = {}
            stat_dict["total"] = float(stat[0])
            stat_dict["hotspot"] = float(stat[1])
            stat_dict["spot"] = float(stat[2])
            spot["stat"][s] = stat_dict
    json_dump(spot)

## convert json files to latex report
def spot_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"samples": "", "reps": ""}):
    data = json_load(input["json"])["stat"]
    spot = []
    for s in param["samples"]:
        spot.append(decimal_to_latex_percent(data[s]["spot"]))

    spot_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "peaks_calling",
        param = {"section_name": "spot",
                 "render_dump": output["latex"],
                 "spot": spot,
                 "reps": param["reps"]})

    template_dump(spot_latex)

## default hotspot peaks calling for all reads and 5M reads replicates files
## SPOT score to json and latex on 5M reads for 2 replicates data
def _hotspot_on_replicates(workflow, conf, tex):
    """
    Use hotspot v4, output list:
	*-final/*.hot.bed               minimally thresholded hotspots( corresponding to hotspot v3 b, broad Peak)
	*-final/*.fdr0.01.hot.bed       FDR thresholded hotspots  ( corresponding to hotspot v3 c)
	*-final/*.fdr0.01.pks.bed       FDR thresholded peaks     ( corresponding to hotspot v3 d, narrow Peak)
	*-final/tag.density.starch      20bp resolution, converted to bigwiggle
    """
    if conf.seq_type.startswith("bed"):
        kind = "_all.bed.starch"
        suffix = "_5M.bed.starch"
    else:
        kind = ".bam"
        suffix = "_5M_sort.bam"

    for i, target in enumerate(conf.treatment_targets):
        ## hotspot v4 support BAM and bed.starch(in a new directory for filtering tags) input
        ## generate configuration for hotspot v4
        hotspot_conf = attach_back(workflow,
            PythonCommand(spot_conf,
                input = {"spot_conf": token_file,
                         ## for bed.starch input files, should be in a different output directory
                         "tag": conf.hotspot_starch_input[i] + kind if conf.seq_type.startswith("bed") else target + kind,
                         "mappable_region": conf.get("hotspot", "mappable_region"),
                         "chrom_info": conf.get("hotspot", "chrom_info")},
                output = {"conf": target + "_runall.tokens.txt",
                          "dir": conf.target_dir},
                param = {"fdrs": "0.01", "K": conf.get("Basis", "read_length"),
                         "species": conf.get("Basis", "species"),
                         "keep_dup": "T",
                         "omit": resource_filename("gcap", "static/Satellite.hg19.bed") if conf.get("Basis", "species") == "hg19" else ""}, ## mouse do not run_badspot
                name = "hotspot config"))
        hotspot_conf.param.update(conf.items("hotspot"))

        ## run hotspot v4
        attach_back(workflow, ShellCommand(
            "{tool} {input[pipe]} {input[token]} {output[dir]} {param[spot]} {param[omit]}",
            tool = "runhotspot",
            input = {"pipe": pipeline_scripts, "token": target + "_runall.tokens.txt"},
            output = {"dir": conf.target_dir,
                      "hotspot": conf.hotspot_reps_final_prefix[i] + ".hot.bed",
                      "density_starch": target + ".tagdensity.bed.starch",
                      "peaks": conf.hotspot_reps_final_prefix[i] + ".fdr0.01.pks.bed"},
            param = {"spot": "all",
                     "omit": resource_filename("gcap", "static/Satellite.hg19.bed") if conf.get("Basis", "species") == "hg19" else ""}))

        ## hotspot all reads bigwiggle
        attach_back(workflow, ShellCommand(
            "{tool} {input[starch]} > {output[bed]}",
            tool = "unstarch",
            input = {"starch": target + ".tagdensity.bed.starch"},
            output = {"bed": target + ".tagdensity.bed.tmp"}))
        attach_back(workflow,
            ShellCommand(
                '{tool} intersect -a {input} -b {param[chrom_bed]} -f 1.00 > {output}',
                tool="bedtools",
                input=target + ".tagdensity.bed.tmp",
                output=target + ".tagdensity.bed",
                param={'chrom_bed': conf.get("lib", "chrom_bed")},
                name="bed replicate filtering"))

        ## convert to bigwiggle, 20 bp resolution
        attach_back(workflow,
            ShellCommand(
                'cut -f 1,2,3,5 {input} > {input}.tmp && {tool} {input}.tmp {param[chrom_len]} {output}',
                tool = "bedGraphToBigWig",
                input=target + ".tagdensity.bed",
                output=target + "_density.bw",
                param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))

    ## Use 5M reads to estimate SPOT
    for i, target in enumerate(conf.treatment_targets):
        ## prepare input for hotspot, get mappable tags and starch into output directory
        if conf.seq_type.startswith("bed"):  ## get library complexity
            attach_back(workflow,
                ShellCommand(
                    "{tool} {input} {output[starch]} {param[tool]} {output[redundancy]}",
                    tool = "bed_duplicates.sh",
                    input = target + "_5M.bed",
                    output = {"starch": conf.hotspot_starch_input[i] + suffix,
                              "redundancy": target + ".dup_metrics"},
                    param = {"tool": conf.peakcalltool}))

        ## generate configuration for hotspot v4
        hotspot_conf = attach_back(workflow,
            PythonCommand(spot_conf,
                input = {"spot_conf": token_file,
                         "tag": conf.hotspot_starch_input[i] + suffix if conf.seq_type.startswith("bed") else target + suffix, ## for bed.starch input files, should be in a different output directory
                         "mappable_region": conf.get("hotspot", "mappable_region"),
                         "chrom_info": conf.get("hotspot", "chrom_info")},
                output = {"conf": target + "_runall_5M.tokens.txt",
                          "dir": conf.target_dir},
                param = {"fdrs": "0.01", "K": conf.get("Basis", "read_length"),
                         "species": conf.get("Basis", "species"),
                         "keep_dup": "T"}))
        hotspot_conf.param.update(conf.items("hotspot"))
        ## run hotspot v4
        hotspot_run = attach_back(workflow, ShellCommand(
            "{tool} {input[pipe]} {input[token]} {output[dir]} {param[spot]} {param[omit]}",
            tool = "runhotspot",
            input = {"pipe": pipeline_scripts, "token": target + "_runall_5M.tokens.txt"},
            output = {"dir": conf.target_dir,
                      "spot_peaks_combined": conf.hotspot_reps_final_5M_prefix[i] + ".fdr0.01.pks.bed",
                      "density_starch": target + "_5M_sort.tagdensity.bed.starch",
                      "hotspot":  conf.hotspot_reps_final_5M_prefix[i] + ".hot.bed",
                      "spot": target + "_5M_sort.spot.out"},
            param = {"spot": "5M",
                     "omit": resource_filename("gcap", "static/Satellite.hg19.bed") if conf.get("Basis", "species") == "hg19" else ""}))

        ## 5M reads density from hotspot v4 tag density, 20bp resolution
        ## for correlation evaluation
        attach_back(workflow, ShellCommand(
            "{tool} {input[starch]} > {output[bed]}",
            tool = "unstarch",
            input = {"starch": target + "_5M_sort.tagdensity.bed.starch"},
            output = {"bed": target + "_5M.tagdensity.bed.tmp"}))
        ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
        ## filter bdg file to remove over-border coordinates
        ## use 5M bigwiggle to estimate replicates consistency
        attach_back(workflow,
            ShellCommand(
                '{tool} intersect -wa -a {input} -b {param[chrom_bed]} -f 1.00 > {output}',
                tool="bedtools",
                input=target + "_5M.tagdensity.bed.tmp",
                output=target + "_5M.tagdensity.bed",
                param={'chrom_bed': conf.get("lib", "chrom_bed")},
                name="bed replicate filtering"))
        ## convert to bigwiggle
        attach_back(workflow,
            ShellCommand(
                'cut -f 1,2,3,5 {input} > {input}.final && {tool} {input}.final {param[chrom_len]} {output}',
                tool = "bedGraphToBigWig",
                input=target + "_5M.tagdensity.bed",
                output=target + "_5M.bw",
                param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))

    ## depend on bed_duplicates.sh, calculate redundancy using awk
    if conf.seq_type.startswith("bed"):
        attach_back(workflow, PythonCommand(
            stat_redun_census,
            input = {"redundancy": [ target + ".dup_metrics" for target in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_redun.json"}, param = {"samples": conf.treatment_bases, "format": conf.seq_type}))
        attach_back(workflow, PythonCommand(
            redundancy_doc,
            input = {"tex": tex, "json": conf.json_prefix + "_redun.json"},
            output = {"redun": conf.latex_prefix + "redun.tex"},
            param = {"reps": len(conf.treatment_pairs), "samples": conf.treatment_bases}))

## call peaks on merged reads files by hotspot
def _hotspot_combo(workflow, conf):
    ## all reads combo
    if conf.seq_type.startswith("bed"):
        merged = conf.hotspot_merge_starch + "_merge_all.bed.starch"
        attach_back(workflow,
            ShellCommand(
                "{tool} -u {param[tag_list]} | sort-bed - | starch - > {output}",
                tool = "bedops",
                input = [ target + ".bed.starch" for target in conf.treatment_targets ],
                output = merged,
                param = {"tag_list": " ".join([ target + ".bed.starch" for target in conf.treatment_targets ])}))
    else:
        merged = conf.prefix + "_merge_all.bam"
        merge_bams_treat = ShellCommand(
            "{tool} merge {output[merged]} {param[bams]}",
            tool="samtools",
            input=[ target + ".bam" for target in conf.treatment_targets],
            output={"merged": merged})
        merge_bams_treat.param = {"bams": " ".join(merge_bams_treat.input)}
        attach_back(workflow, merge_bams_treat)

    ## config for hotspot v4
    hotspot_conf = attach_back(workflow,
        PythonCommand(spot_conf,
            input = {"spot_conf": token_file,
                     "tag": merged,
                     "mappable_region": conf.get("hotspot", "mappable_region"),
                     "chrom_info": conf.get("hotspot", "chrom_info")},
            output = {"conf": conf.prefix + "_runall.tokens.txt",
                      "dir": conf.target_dir},
            param = {"fdrs": "0.01", "K": conf.get("Basis", "read_length"),
                     "species": conf.get("Basis", "species"),
                     "keep_dup": "T",
                     "omit": resource_filename("gcap", "static/Satellite.hg19.bed") if conf.get("Basis", "species") == "hg19" else ""}))
    hotspot_conf.param.update(conf.items("hotspot"))

    ## run hotspot v3
    attach_back(workflow, ShellCommand(
        "{tool} {input[pipe]} {input[token]} {output[dir]} {param[spot]} {param[omit]}",
        tool = "runhotspot",
        input = {"pipe": pipeline_scripts, "token": conf.prefix + "_runall.tokens.txt"},
        output = {"dir": conf.target_dir,
                  "hotspot": conf.hotspot_merge_final_prefix + "_merge_all.hot.bed",
                  "peaks": conf.hotspot_merge_final_prefix + "_merge_all.fdr0.01.pks.bed",
                  "density_starch": conf.prefix + "_merge_all.tagdensity.bed.starch"},
        param = {"spot": "all",
                 "omit": resource_filename("gcap", "static/Satellite.hg19.bed") if conf.get("Basis", "species") == "hg19" else ""}))
    attach_back(workflow, ShellCommand(
        "{tool} {input[starch]} > {output[bed]}",
        tool = "unstarch",
        input = {"starch": conf.prefix + "_merge_all.tagdensity.bed.starch"},
        output = {"bed": conf.prefix + "_merge_all.tagdensity.bed.tmp"},
        param={'chrom_bed': conf.get("lib", "chrom_bed")},
        name="bed replicate filtering"))
    attach_back(workflow,
        ShellCommand(
            '{tool} intersect -a {input} -b {param[chrom_bed]} -f 1.00 > {output}',
            tool="bedtools",
            input=conf.prefix + "_merge_all.tagdensity.bed.tmp",
            output=conf.prefix + "_merge_all.tagdensity.bed",
            param={'chrom_bed': conf.get("lib", "chrom_bed")},
            name="bed replicate filtering"))
    ## convert to bigwiggle
    attach_back(workflow,
        ShellCommand(
            'cut -f 1,2,3,5 {input} > {input}.final && {tool} {input}.final {param[chrom_len]} {output}',
            tool = "bedGraphToBigWig",
            input=conf.prefix + "_merge_all.tagdensity.bed",
            output=conf.prefix + ".bw",
            param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))

## call peaks on replicates reads files by macs2
def _macs2_on_reps(workflow, conf, tex):
    """
    call peaks by MACS2, optionally
    """
    ## for all reads
    if conf.seq_type.startswith("bed"):
        kind = "_all.bed"
        suffix = "_5M.bed"
    else:
        kind = ".bam"
        suffix = "_5M_sort.bam"

    for target in conf.treatment_targets:
        ## keep all duplicate tags as hotspot
        macs2_on_rep = attach_back(workflow,
            ShellCommand(
                "{tool} callpeak -B -q {param[fdr]} -f {param[format]} -g {param[species]} --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel -g {param[species]} \
                {param[treat_opt]} -n {param[description]}",
                tool="macs2",
                input={"treat": target + kind},
                output={"peaks": target + "_macs2_all_peaks.bed",
                        "summit": target + "_macs2_all_summits.bed",
                        "treat_bdg": target + "_macs2_all_treat_pileup.bdg",
                        "peaks_xls": target + "_macs2_all_peaks.xls",
                        "control_bdg": target + "_macs2_all_control_lambda.bdg"},
                param={"description": target + "_macs2_all", "keep_dup": "all", "shiftsize": 50, "species": conf.get("macs2", "species"),
                       "fdr": 0.01},
                name="macs2_callpeak_rep"))
        macs2_on_rep.param["treat_opt"] = "-t " + macs2_on_rep.input["treat"]
        macs2_on_rep.update(param=conf.items("macs2"))
        if conf.seq_type == "se":
            macs2_on_rep.param["format"] = "BAM"
        elif conf.seq_type == "pe":
            macs2_on_rep.param["format"] = "BAMPE"   ## pair end mode, not compatible with SAM files built-in sampling
        elif conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam"):
            if conf.seq_type.split(",")[1].strip().lower() == "se":
                macs2_on_rep.param["format"] = "BAM"
            elif conf.seq_type.split(",")[1].strip().lower() == "pe":
                macs2_on_rep.param["format"] = "BAMPE"
        elif conf.seq_type.startswith("bed"):
            macs2_on_rep.param["format"] = "BED"

        ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
        ## filter bdg file to remove over-border coordinates
        attach_back(workflow,
            ShellCommand(
                '{tool} intersect -a {input} -b {param[chrom_bed]} -f 1.00 > {output}',
                tool="bedtools",
                input=target + "_macs2_all_treat_pileup.bdg",
                output=target + "_macs2_all_treat_pileup.bdg.tmp",
                param={'chrom_bed': conf.get("lib", "chrom_bed")},
                name="bedGraph replicate filtering"))
        bdg2bw_treatrep = attach_back(workflow,
            ShellCommand(
                "{tool} {input} {param[chrom_len]} {output}",
                tool="bedGraphToBigWig",
                input=target + "_macs2_all_treat_pileup.bdg.tmp",
                output=target + "_treat_all.bw",
                param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))

    ## for 5M reads
    if conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam"):
        ty = conf.seq_type.split(",")[1].strip().lower()
    else:
        ty = conf.seq_type

    for target in conf.treatment_targets:
        ## keep all duplicate tags as hotspot
        macs2_on_rep = attach_back(workflow,
            ShellCommand(
                "{tool} callpeak -B -q {param[fdr]} -g {param[species]} --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel -g {param[species]} \
                {param[treat_opt]} -n {param[description]}",
                tool="macs2",
                input={"treat": target + suffix},
                output={"peaks": target + "_5M_macs2_peaks.bed",
                        "summit": target + "_5M_macs2_summits.bed",
                        "treat_bdg": target + "_5M_macs2_treat_pileup.bdg",
                        "peaks_xls": target + "_5M_macs2_peaks.xls",
                        "control_bdg": target + "_5M_macs2_control_lambda.bdg"},
                param={"description": target + "_5M_macs2", "keep_dup": "all", "shiftsize": 50, "species": conf.get("macs2", "species"),
                       "fdr": 0.01},
                name="macs2_callpeak_rep"))
        macs2_on_rep.param["treat_opt"] = "-t " + macs2_on_rep.input["treat"]
        macs2_on_rep.update(param=conf.items("macs2"))

        ## redundancy ratio for reads BED files
        if conf.seq_type.startswith("bed"):
            attach_back(workflow,
                ShellCommand(
                    "{tool} {input} {param[starch]} {param[tool]} {output[redundancy]}",
                    tool = "bed_duplicates.sh",
                    input = target + suffix,
                    output = {"redundancy": target + ".dup_metrics"},
                    param = {"tool": conf.peakcalltool,
                             "starch": "no"})) ## no need to starch for macs2

        ## SPOT score for MACS2 5M reads
        attach_back(workflow, ShellCommand(
            "{tool} {input[starch]} {input[bed]}",
            tool = "macs2_spot.sh",
            input = {"starch": target + suffix, "bed": target + "_5M_macs2_peaks.bed"},
            output = target + "_5M_macs2_peaks.bed" + ".spot.out"))

        ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
        ## filter bdg file to remove over-border coordinates
        attach_back(workflow,
            ShellCommand(
                '{tool} intersect -a {input} -b {param[chrom_bed]} -f 1.00 > {output}',
                tool="bedtools",
                input=target + "_5M_macs2_treat_pileup.bdg",
                output=target + "_5M_macs2_treat_pileup.bdg.tmp",
                param={'chrom_bed': conf.get("lib", "chrom_bed")},
                name="bedGraph replicate filtering"))
        bdg2bw_treatrep = attach_back(workflow,
            ShellCommand(
                "{tool} {input} {param[chrom_len]} {output}",
                tool="bedGraphToBigWig",
                input=target + "_5M_macs2_treat_pileup.bdg.tmp",
                output=target + "_5M_macs2_treat.bw",
                param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))

    if conf.seq_type.startswith("bed"):
        attach_back(workflow, PythonCommand(
            stat_redun_census,
            input = {"redundancy": [ target + ".dup_metrics" for target in conf.treatment_targets ]},
            output = {"json": conf.json_prefix + "_redun.json"}, param = {"samples": conf.treatment_bases, "format": conf.seq_type}))
        attach_back(workflow, PythonCommand(
            redundancy_doc,
            input = {"tex": tex, "json": conf.json_prefix + "_redun.json"},
            output = {"redun": conf.latex_prefix + "redun.tex"},
            param = {"reps": len(conf.treatment_pairs), "samples": conf.treatment_bases}))

## call peaks by macs2 on merged reads files
def _macs2_on_combo(workflow, conf):
    # merge all treatments into one
    if not conf.seq_type.startswith("bed"):
        merge_bams_treat = ShellCommand(
            "{tool} merge {output[merged]} {param[bams]}",
            tool="samtools",
            input=[ target + ".bam" for target in conf.treatment_targets],
            output={"merged": conf.prefix + "_treatment.bam"})
        merge_bams_treat.param = {"bams": " ".join(merge_bams_treat.input)}
        attach_back(workflow, merge_bams_treat)
        kind = "_treatment.bam"
    else:
        merge_beds_treat = ShellCommand(
            "{tool} {param[bed]} > {output[bed]}",
            tool = "cat",
            input = [ target + ".bed.all" for target in conf.treatment_targets ],
            output = {"bed": conf.prefix + "_treatment.bed"})
        merge_beds_treat.param = {"bed": " ".join(merge_beds_treat.input)}
        attach_back(workflow, merge_beds_treat)
        kind = "_treatment.bed"
    macs2_on_merged = attach_back(workflow, ShellCommand(
        "{tool} callpeak -B -q {param[fdr]} -f {param[format]} --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel -g {param[species]} \
        {param[treat_opt]} -n {param[description]}",
        tool="macs2",
        input={"merged": conf.prefix + kind},
        output={"peaks": conf.prefix + "_peaks.bed",
                "summit": conf.prefix + "_summits.bed",
                "treat_bdg": conf.prefix + "_treat_pileup.bdg",
                "peaks_xls": conf.prefix + "_peaks.xls",
                "control_bdg": conf.prefix + "_control_lambda.bdg"},
        param={"description": conf.prefix,
               "keep_dup": "all",
               "shiftsize": 50,
               "fdr": 0.01,
               "species": conf.get("macs2",  "species")},
        name="macs2_callpeak_merged"))

    macs2_on_merged.param["treat_opt"] = " -t " + macs2_on_merged.input["merged"]
    macs2_on_merged.update(param=conf.items("macs2"))

    if conf.seq_type == "se":
        macs2_on_merged.param["format"] = "BAM"
    elif conf.seq_type == "pe":
        macs2_on_merged.param["format"] = "BAMPE"   ## pair end mode, not compatible with SAM files built-in sampling
    elif conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam"):
        if conf.seq_type.split(",")[1].strip().lower() == "se":
            macs2_on_merged.param["format"] = "BAM"
        elif conf.seq_type.split(",")[1].strip().lower() == "pe":
            macs2_on_merged.param["format"] = "BAMPE"
    elif conf.seq_type.startswith("bed"):
        macs2_on_merged.param["format"] = "BED"

    ## BedClip
    ## prototype used here to do the similar thing on bedclip
    bed_clip = attach_back(workflow,
        ShellCommand(
            template="{tool} {input} {param[chrom_len]} {output}",
            tool="bedClip",
            input=conf.prefix + "_treat_pileup.bdg",
            output =conf.prefix + "_treat_pileup.bdg.tmp",
            param={'chrom_len': conf.get_path("lib", "chrom_len")},
            name="bedclip filter"))

    bdg2bw_treat = attach_back(workflow,
        ShellCommand(
            "{tool} {input[bdg]} {input[chrom_len]} {output[bw]}",
            tool="bedGraphToBigWig",
            input={"bdg": conf.prefix + "_treat_pileup.bdg.tmp",
                   "chrom_len": conf.get("lib", "chrom_len")},
            output={"bw": conf.prefix + "_macs2_treat_all.bw"},
            name="bdg_to_bw"))

def _peaks_calling_latex(workflow, conf, tex):
    ## peaks number json
    if conf.peakcalltool == "hotspot":
        if len(conf.treatment_pairs) >= 2:
            ## use hotspot b to evaluate 5M reads peaks number, replicates consistency
            ## use peaks d as narrow peaks to evaluate all reads peaks number
            attach_back(workflow, PythonCommand(
                stat_peaks,
                input = {"peaks": {"5M_spot": [ conf.hotspot_reps_final_5M_prefix[i] + ".hot.bed" for i in range(len(conf.treatment_pairs)) ],
                                   "all_peaks": [ conf.hotspot_reps_final_prefix[i] + ".fdr0.01.pks.bed" for i in range(len(conf.treatment_pairs)) ],
                                   "combo": conf.hotspot_merge_final_prefix + "_merge_all.fdr0.01.pks.bed" }},
                output = {"json": conf.json_prefix + "_peaks.json"},
                param = {"tool": conf.peakcalltool, "samples": conf.treatment_bases},
                name = "json peaks and hotspot"))
        else:
            attach_back(workflow, PythonCommand(
                stat_peaks,
                input = {"peaks": {"5M_spot": [ conf.hotspot_reps_final_5M_prefix[i] + ".hot.bed" for i in range(len(conf.treatment_pairs)) ],
                                   "all_peaks": [ conf.hotspot_reps_final_prefix[i] + ".fdr0.01.pks.bed" for i in range(len(conf.treatment_pairs)) ]}},
                output = {"json": conf.json_prefix + "_peaks.json"},
                param = {"tool": conf.peakcalltool, "samples": conf.treatment_bases},
                name = "json peaks and hotspot"))

    elif conf.peakcalltool == "macs2":
        if len(conf.treatment_pairs) >= 2:
            attach_back(workflow, PythonCommand(
                stat_peaks,
                input = {"peaks": {"all_peaks": [ target + "_macs2_all_peaks.bed" for target in conf.treatment_targets ],
                                   "5M_spot": [ target + "_5M_macs2_peaks.bed" for target in conf.treatment_targets ],
                                   "combo": conf.prefix + "_peaks.bed"}},
                output = {"json": conf.json_prefix + "_peaks.json"},
                param = {"tool": conf.peakcalltool, "samples": conf.treatment_bases}))
        else:
            attach_back(workflow, PythonCommand(
                stat_peaks,
                input = {"peaks": {"all_peaks": [ target + "_macs2_all_peaks.bed" for target in conf.treatment_targets ],
                                   "5M_spot": [ target + "_5M_macs2_peaks.bed" for target in conf.treatment_targets ]}},
                output = {"json": conf.json_prefix + "_peaks.json"},
                param = {"tool": conf.peakcalltool, "samples": conf.treatment_bases}))

    ## peaks number latex
    attach_back(workflow, PythonCommand(
        peaks_doc,
        input = {"tex": tex, "json": conf.json_prefix + "_peaks.json"},
        output = {"latex": conf.latex_prefix + "_peaks.tex"},
        param = {"reps": len(conf.treatment_pairs), "samples": conf.treatment_bases},
        name = "peaks number json"))

    if conf.peakcalltool == "hotspot":
        spot = [ target + "_5M_sort.spot.out" for target in conf.treatment_targets ]
    elif conf.peakcalltool == "macs2":
        spot = [ target + "_5M_macs2_peaks.bed" + ".spot.out" for target in conf.treatment_targets ]

    ## 5M reads, calculate the merged two passes and peaks regions number
    ## 5M reads for macs2 optionally
    attach_back(workflow, PythonCommand(
        stat_spot_on_replicates,
        input = {"spot_files": spot},
        output = {"json": conf.json_prefix + "_sample_spot_5M.json"},
        name = "json spot",
        param = {"samples": conf.treatment_bases}))

    attach_back(workflow, PythonCommand(
        spot_doc,
        input = {"json": conf.json_prefix + "_sample_spot_5M.json",
                 "tex": tex},
        output = {"latex": conf.latex_prefix + "_spot.tex"},
        name = "doc spot 5M",
        param = {"samples": conf.treatment_bases, "reps": len(conf.treatment_pairs)}))

def stat_peaks(input = {"peaks": {"all_peaks": "", "5M_spot": "", "combo_peaks": ""}}, output = {"json"}, param = {"tool": "", "samples": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}

    json_dict["stat"]["all_peaks"] = {}

    # d, narrow peaks
    for a, s in zip(input["peaks"]["all_peaks"], param["samples"]):
        with open(a) as f:
            json_dict["stat"]["all_peaks"][s] = len(f.readlines())

    json_dict["stat"]["5M_spot"] = {}

    # b, hotspot minimally threshold
    for spot, sam in zip(input["peaks"]["5M_spot"], param["samples"]):
        with open(spot) as f:
            json_dict["stat"]["5M_spot"][sam] = len(f.readlines()[1:])

    json_dict["stat"]["combo"] = {}

    if len(input["peaks"]["all_peaks"]) >= 2:
        with open(input["peaks"]["combo"]) as combo:
            json_dict["stat"]["combo"]= len(combo.readlines())
    print(json_dict)

    json_dump(json_dict)

def peaks_doc(input = {"tex": "", "json": ""}, output = {"latex": ""}, param = {"samples": "", "reps": ""}):
    """ use d peaks( narrow peaks ) for all reads evaluation
     use b Hotspot for 5M reads evaluation and replicates consistency(correlation and overlap)
     use 5M narrow peaks for DHS evaluation
    """
    data = json_load(input["json"])
    peaks_5M = []
    peaks_all = []

    for s in param["samples"]:
        peaks_5M.append(data["stat"]["5M_spot"][s])

    for a in param["samples"]:
        peaks_all.append(data["stat"]["all_peaks"][a])

    if data["stat"]["combo"]:
        combo = data["stat"]["combo"]
        peaks_latex = JinjaTemplateCommand(
            template = input["tex"],
            name = "peaks_calling",
            param = {"section_name": "peaks_calling",
                     "render_dump": output["latex"],
                     "spot_5M": peaks_5M,
                     "peaks_all":peaks_all, "combo": combo,
                     "tool": data["param"]["tool"],
                     "reps": param["reps"]})

    else:
        peaks_latex = JinjaTemplateCommand(
            template = input["tex"],
            name = "peaks_calling",
            param = {"section_name": "peaks_calling",
                     "render_dump": output["latex"],
                     "spot_5M": peaks_5M,
                     "peaks_all": peaks_all,
                     "tool": data["param"]["tool"],
                     "reps": param["reps"]})
    template_dump(peaks_latex)