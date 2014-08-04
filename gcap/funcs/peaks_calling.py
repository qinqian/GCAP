
###############################################################
#
# 1. hotspot v4 peak calling
# 5. peaks number evaluation
# 6. replicates consistency by UCSC utility, edwReplicatePeak, bigwigCorrelate
#
###############################################################

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back
from gcap.funcs.helpers import *
from pkg_resources import resource_filename

## call peaks by hotspot or macs2 to QC DNase data
def call_peaks(workflow, conf, tex):

    if 'hotspot' in conf.get("tool", "peak_calling"):
        hotspotv4(workflow, conf, tex)
        # if conf.peakcalltool == "hotspot":

        # elif conf.peakcalltool == "macs2":
        #     _macs2_on_reps(workflow, conf, tex)
        ## evaluate replicates consistency
        # peaks_reps_evaluating(workflow, conf, tex)
        ## peaks calling on merged reads files
        # if conf.peakcalltool == "hotspot":
        #     _hotspot_combo(workflow, conf)
        # elif conf.peakcalltool == "macs2":
        #     _macs2_on_combo(workflow, conf)

def eval_reps(workflow, conf, tex):

    attach_back(workflow, ShellCommand(
        """
        cat {param[narrowPeaks]} | sort -k1,1 -k2,2n - | bedtools merge -i - > {output[mergedPeak]}
        bedToBigBed {output[mergedPeak]} {param[chromsize]} {output[mergedPeakbb]}
        bigWigCorrelate -restrict={output[mergedPeakbb]} {param[bigwigs]} 1>{output[qc1]}
        {tool} {param[narrowPeaksbb]} {output[qc2]}
        """,
        tool = "edwComparePeaks",
        input = {"narrowPeaks": [ target + ".narrowPeak" for target in conf.treatment_targets ],
                 "bigwigs": [ target + ".bigWig" for target in conf.treatment_targets ],
                 "narrowPeakbbs": [ target + ".narrowPeak.bigBed" for target in conf.treatment_targets ]},
        output = {"mergedPeak": conf.prefix + "_merge.bed",
                  "mergedPeakbb": conf.prefix + "_merged.bigBed",
                  "qc1": conf.prefix + "_cor.qc",
                  "qc2": conf.prefix + "_overlap.qc"},
        param = {"narrowPeaksbb": " ".join([ target + ".narrowPeak.bigBed" for target in conf.treatment_targets ]),
                 "narrowPeaks": " ".join([ target + ".narrowPeak" for target in conf.treatment_targets ]),
                 "bigwigs": " ".join([ target + ".bigWig" for target in conf.treatment_targets ]),
                 "chromsize": conf.get(conf.species, "chrom_len")}))


def hotspotv4(workflow, conf, tex):
    for target in conf.treatment_targets:
        narrow = resource_filename("gcap", "glue/narrowPeak.as")
        broad = resource_filename("gcap", "glue/broadPeak.as")
        hotspot=attach_back(workflow,
                    ShellCommand(
                        "{tool} {param[hotspot_dir]} {param[genome]} {input[bam]} {param[readsize]} {output[narrowbb]} {output[broadbb]} {output[bigwig]} {param[tmp]} {output[hotspot_output]} {input[narrowas]} {input[broadas]} {param[chromsize]} {output[narrow]} {output[broad]}",
                        tool = "eap_run_hotspot",
                        input = {"bam": target + "_final_nochrm.bam",
                                 "narrowas": narrow,
                                 "broadas": broad},
                        output = {"narrowbb": target + ".narrowPeak.bigBed",
                                  "broadbb": target + ".broadPeak.bigBed",
                                  "narrow": target + ".narrowPeak",
                                  "broad": target + ".broadPeak",
                                  "bigwig": target + ".bigWig",
                                  "hotspot_output": target + "_hotspot"},
                        param = {"hotspot_dir": conf.get("tool", "peak_calling"),
                                 "genome": conf.species,
                                 "chromsize": conf.get(conf.species, "chrom_len"),
                                 "tmp": target + "_hotspot_peak_call_tmp",
                                 "readsize": 36}))
    have_treat_reps = len(conf.treatment_pairs) >= 2 ## replicates

    if have_treat_reps:
        eval_reps(workflow, conf, tex)
        catsam = attach_back(workflow, ShellCommand(
            "{tool} cat {param[bams]} > {output[bam]}",
            tool = "samtools",
            input ={"bams": [ target + "_final.bam" for target in conf.treatment_targets]},
            output = {"bam": conf.prefix + "_pool.bam"}))
        catsam.param.update(bams=' '.join(catsam.input["bams"]))
        hotspot_merge = hotspot.clone
        hotspot_merge.param.update(tmp=conf.prefix+"_hotspot_peak_call_tmp")
        hotspot_merge.input.update(bam = conf.prefix + "_pool.bam")
        hotspot_merge.output ={"narrowbb": conf.prefix + ".narrowPeak.bigBed",
                               "broadbb": conf.prefix + ".broadPeak.bigBed",
                               "narrow": conf.prefix + ".narrowPeak",
                               "broad": conf.prefix + ".broadPeak",
                               "bigwig": conf.prefix + ".bigWig",
                               "hotspot_output": conf.prefix + "_hotspot"}
        attach_back(workflow, hotspot_merge)

# ## default hotspot peaks calling for all reads and 5M reads replicates files
# ## SPOT score to json and latex on 5M reads for 2 replicates data
# def _hotspot_on_replicates(workflow, conf, tex):
#     """
#     """
#     if conf.seq_type.startswith("bed"):
#         kind = ".bed.starch"
#         suffix = "_5M_sort.bed.starch"
#     else:
#         kind = ".bam"
#         suffix = "_5M_sort.bam"

#     for i, target in enumerate(conf.treatment_targets):
#         ## hotspot v4 support BAM and bed.starch(in a new directory for filtering tags) input
#         ## generate configuration for hotspot v4
#         hotspot_conf = attach_back(workflow,
#             PythonCommand(spot_conf,
#                 input = {"spot_conf": token_file,
#                          ## for bed.starch input files, should be in a different output directory
#                          "tag": conf.hotspot_starch_input[i] + kind if conf.seq_type.startswith("bed") else target + kind,
#                          "mappable_region": conf.get("hotspot", "mappable_region"),
#                          "chrom_info": conf.get("hotspot", "chrom_info")},
#                 output = {"conf": target + "_runall.tokens.txt",
#                           "dir": conf.target_dir},
#                 param = {"fdrs": "0.01", "K": conf.get("Basis", "read_length"),
#                          "species": conf.get("Basis", "species"),
#                          "keep_dup": "T",
#                          "omit": resource_filename("gcap", "static/Satellite.hg19.bed") if conf.get("Basis", "species") == "hg19" else ""}, ## mouse do not run_badspot
#                 name = "hotspot config"))
#         hotspot_conf.param.update(conf.items("hotspot"))

# ## call peaks on merged reads files by hotspot
# def _hotspot_combo(workflow, conf):
#     ## all reads combo
#     if conf.seq_type.startswith("bed"):
#         merged = conf.hotspot_merge_starch + "_merge_all.bed.starch"
#         attach_back(workflow,
#             ShellCommand(
#                 "{tool} -u {param[tag_list]} | sort-bed - | starch - > {output}",
#                 tool = "bedops",
#                 input = [ target + ".bed.starch" for target in conf.treatment_targets ],
#                 output = merged,
#                 param = {"tag_list": " ".join([ target + ".bed.starch" for target in conf.treatment_targets ])}))
#     else:
#         merged = conf.prefix + "_merge_all.bam"
#         merge_bams_treat = ShellCommand(
#             "{tool} merge {output[merged]} {param[bams]}",
#             tool="samtools",
#             input=[ target + ".bam" for target in conf.treatment_targets],
#             output={"merged": merged})
#         merge_bams_treat.param = {"bams": " ".join(merge_bams_treat.input)}
#         attach_back(workflow, merge_bams_treat)

#     ## config for hotspot v4
#     hotspot_conf = attach_back(workflow,
#         PythonCommand(spot_conf,
#             input = {"spot_conf": token_file,
#                      "tag": merged,
#                      "mappable_region": conf.get("hotspot", "mappable_region"),
#                      "chrom_info": conf.get("hotspot", "chrom_info")},
#             output = {"conf": conf.prefix + "_runall.tokens.txt",
#                       "dir": conf.target_dir},
#             param = {"fdrs": "0.01", "K": conf.get("Basis", "read_length"),
#                      "species": conf.get("Basis", "species"),
#                      "keep_dup": "T",
#                      "omit": resource_filename("gcap", "static/Satellite.hg19.bed") if conf.get("Basis", "species") == "hg19" else ""}))
#     hotspot_conf.param.update(conf.items("hotspot"))

#     ## run hotspot v3
#     attach_back(workflow, ShellCommand(
#         "{tool} {input[pipe]} {input[token]} {output[dir]} {param[spot]} {param[omit]}",
#         tool = "runhotspot",
#         input = {"pipe": pipeline_scripts, "token": conf.prefix + "_runall.tokens.txt"},
#         output = {"dir": conf.target_dir,
#                   "hotspot": conf.hotspot_merge_final_prefix + "_merge_all.hot.bed",
#                   "peaks": conf.hotspot_merge_final_prefix + "_merge_all.fdr0.01.pks.bed",
#                   "density_starch": conf.prefix + "_merge_all.tagdensity.bed.starch"},
#         param = {"spot": "all",
#                  "omit": resource_filename("gcap", "static/Satellite.hg19.bed") if conf.get("Basis", "species") == "hg19" else ""}))
#     attach_back(workflow, ShellCommand(
#         "{tool} {input[starch]} > {output[bed]}",
#         tool = "unstarch",
#         input = {"starch": conf.prefix + "_merge_all.tagdensity.bed.starch"},
#         output = {"bed": conf.prefix + "_merge_all.tagdensity.bed.tmp"},
#         param={'chrom_bed': conf.get("lib", "chrom_bed")},
#         name="bed replicate filtering"))
#     attach_back(workflow,
#         ShellCommand(
#             '{tool} intersect -a {input} -b {param[chrom_bed]} -f 1.00 > {output}',
#             tool="bedtools",
#             input=conf.prefix + "_merge_all.tagdensity.bed.tmp",
#             output=conf.prefix + "_merge_all.tagdensity.bed",
#             param={'chrom_bed': conf.get("lib", "chrom_bed")},
#             name="bed replicate filtering"))
#     ## convert to bigwiggle
#     attach_back(workflow,
#         ShellCommand(
#             'cut -f 1,2,3,5 {input} > {input}.final && {tool} {input}.final {param[chrom_len]} {output}',
#             tool = "bedGraphToBigWig",
#             input=conf.prefix + "_merge_all.tagdensity.bed",
#             output=conf.prefix + ".bw",
#             param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))

# ## call peaks on replicates reads files by macs2
# def _macs2_on_reps(workflow, conf, tex):
#     """
#     call peaks by MACS2, optionally
#     """
#     ## for all reads
#     if conf.seq_type.startswith("bed"):
#         kind = "_all.bed"
#         suffix = "_5M.bed"
#     else:
#         kind = ".bam"
#         suffix = "_5M_sort.bam"

#     for target in conf.treatment_targets:
#         ## keep all duplicate tags as hotspot
#         macs2_on_rep = attach_back(workflow,
#             ShellCommand(
#                 "{tool} callpeak -B -q {param[fdr]} -f {param[format]} --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel -g {param[species]} \
#                 {param[treat_opt]} -n {param[description]}",
#                 tool="macs2",
#                 input={"treat": target + kind},
#                 output={"peaks": target + "_macs2_all_peaks.encodePeak",
#                         "summit": target + "_macs2_all_summits.bed",
#                         "treat_bdg": target + "_macs2_all_treat_pileup.bdg",
#                         "peaks_xls": target + "_macs2_all_peaks.xls",
#                         "control_bdg": target + "_macs2_all_control_lambda.bdg"},
#                 param={"description": target + "_macs2_all", "keep_dup": "all", "shiftsize": 50, "species": conf.get("macs2", "species"),
#                        "fdr": 0.01},
#                 name="macs2_callpeak_rep"))
#         macs2_on_rep.param["treat_opt"] = "-t " + macs2_on_rep.input["treat"]
#         macs2_on_rep.update(param=conf.items("macs2"))
#         if conf.seq_type == "se":
#             macs2_on_rep.param["format"] = "BAM"
#         elif conf.seq_type == "pe":
#             macs2_on_rep.param["format"] = "BAMPE"   ## pair end mode, not compatible with SAM files built-in sampling
#         elif conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam"):
#             if conf.seq_type.split(",")[1].strip().lower() == "se":
#                 macs2_on_rep.param["format"] = "BAM"
#             elif conf.seq_type.split(",")[1].strip().lower() == "pe":
#                 macs2_on_rep.param["format"] = "BAMPE"
#         elif conf.seq_type.startswith("bed"):
#             macs2_on_rep.param["format"] = "BED"

#         ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
#         ## filter bdg file to remove over-border coordinates
#         attach_back(workflow,
#             ShellCommand(
#                 '{tool} intersect -a {input} -b {param[chrom_bed]} -f 1.00 > {output}',
#                 tool="bedtools",
#                 input=target + "_macs2_all_treat_pileup.bdg",
#                 output=target + "_macs2_all_treat_pileup.bdg.tmp",
#                 param={'chrom_bed': conf.get("lib", "chrom_bed")},
#                 name="bedGraph replicate filtering"))
#         bdg2bw_treatrep = attach_back(workflow,
#             ShellCommand(
#                 "{tool} {input} {param[chrom_len]} {output}",
#                 tool="bedGraphToBigWig",
#                 input=target + "_macs2_all_treat_pileup.bdg.tmp",
#                 output=target + "_treat_all.bw",
#                 param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))

#     ## for 5M reads
#     if conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam"):
#         ty = conf.seq_type.split(",")[1].strip().lower()
#     else:
#         ty = conf.seq_type

#     for target in conf.treatment_targets:
#         ## keep all duplicate tags as hotspot
#         macs2_on_rep = attach_back(workflow,
#             ShellCommand(
#                 "{tool} callpeak -B -q {param[fdr]} --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel -g {param[species]} \
#                 {param[treat_opt]} -n {param[description]}",
#                 tool="macs2",
#                 input={"treat": target + suffix},
#                 output={"peaks": target + "_5M_macs2_peaks.encodePeak",
#                         "summit": target + "_5M_macs2_summits.bed",
#                         "treat_bdg": target + "_5M_macs2_treat_pileup.bdg",
#                         "peaks_xls": target + "_5M_macs2_peaks.xls",
#                         "control_bdg": target + "_5M_macs2_control_lambda.bdg"},
#                 param={"description": target + "_5M_macs2", "keep_dup": "all", "shiftsize": 50, "species": conf.get("macs2", "species"),
#                        "fdr": 0.01},
#                 name="macs2_callpeak_rep"))
#         macs2_on_rep.param["treat_opt"] = "-t " + macs2_on_rep.input["treat"]
#         macs2_on_rep.update(param=conf.items("macs2"))

#         ## redundancy ratio for reads BED files
#         if conf.seq_type.startswith("bed"):
#             attach_back(workflow,
#                 ShellCommand(
#                     "{tool} {input} {param[starch]} {param[tool]} {output[redundancy]}",
#                     tool = "bed_duplicates.sh",
#                     input = target + suffix,
#                     output = {"redundancy": target + ".dup_metrics"},
#                     param = {"tool": conf.peakcalltool,
#                              "starch": "no"})) ## no need to starch for macs2

#         ## SPOT score for MACS2 5M reads
#         attach_back(workflow, ShellCommand(
#             "{tool} {input[starch]} {input[bed]}",
#             tool = "dac_macs2_spot.sh",
#             input = {"starch": target + suffix, "bed": target + "_5M_macs2_peaks.encodePeak"},
#             output = target + "_5M_macs2_peaks.encodePeak" + ".spot.out"))

#         ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
#         ## filter bdg file to remove over-border coordinates
#         attach_back(workflow,
#             ShellCommand(
#                 '{tool} intersect -a {input} -b {param[chrom_bed]} -f 1.00 > {output}',
#                 tool="bedtools",
#                 input=target + "_5M_macs2_treat_pileup.bdg",
#                 output=target + "_5M_macs2_treat_pileup.bdg.tmp",
#                 param={'chrom_bed': conf.get("lib", "chrom_bed")},
#                 name="bedGraph replicate filtering"))
#         bdg2bw_treatrep = attach_back(workflow,
#             ShellCommand(
#                 "{tool} {input} {param[chrom_len]} {output}",
#                 tool="bedGraphToBigWig",
#                 input=target + "_5M_macs2_treat_pileup.bdg.tmp",
#                 output=target + "_5M_macs2_treat.bw",
#                 param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))

#     if conf.seq_type.startswith("bed"):
#         attach_back(workflow, PythonCommand(
#             stat_redun_census,
#             input = {"redundancy": [ target + ".dup_metrics" for target in conf.treatment_targets ]},
#             output = {"json": conf.json_prefix + "_redun.json"}, param = {"samples": conf.treatment_bases, "format": conf.seq_type}))
#         attach_back(workflow, PythonCommand(
#             redundancy_doc,
#             input = {"tex": tex, "json": conf.json_prefix + "_redun.json"},
#             output = {"redun": conf.latex_prefix + "redun.tex"},
#             param = {"reps": len(conf.treatment_pairs), "samples": conf.treatment_bases}))

