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
    if "macs" in conf.get("tool", "peak_calling"):
        macs2(workflow, conf, tex)

narrow = resource_filename("gcap", "glue/narrowPeak.as")
broad = resource_filename("gcap", "glue/broadPeak.as")

def eval_reps(workflow, conf, tex):
    peaks = [ target + ".narrowPeak" for target in conf.treatment_targets ]

    attach_back(workflow, ShellCommand(
        """
        cat {param[narrowPeaks]} | sort -k1,1 -k2,2n - | bedtools merge -i - > {output[mergedPeak]}
        bedToBigBed {output[mergedPeak]} {param[chromsize]} {output[mergedPeakbb]}
        bigWigCorrelate -restrict={output[mergedPeakbb]} {param[bigwigs]} 1>{output[qc1]}
        {tool} {param[narrowPeaksbb]} {output[qc2]}
        """,
        tool = "edwComparePeaks",
        input = {"narrowPeaks": peaks,
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
                                  # "qc1": target + ".narrowPeak.qc",
                                  # "qc2": target + ".broadPeak.qc",
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
                               # "qc1": conf.prefix + ".narrowPeak.qc",
                               # "qc2": conf.prefix + ".broadPeak.qc",
                               "broad": conf.prefix + ".broadPeak",
                               "bigwig": conf.prefix + ".bigWig",
                               "hotspot_output": conf.prefix + "_hotspot"}
        attach_back(workflow, hotspot_merge)


# ## call peaks on replicates reads files by macs2
def macs2(workflow, conf, tex):
    """
    call peaks by MACS2, optionally
    """
    for target in conf.treatment_targets:
        ## note: trim 5 column signal value to 1000
        macs2_on_rep = attach_back(workflow,
            ShellCommand(
                """
                if [ ! -s {output[narrow]} ];then
                {param[tool]} callpeak --SPMR -B -q {param[fdr]} {param[format]} --keep-dup {param[keep_dup]} --extsize={param[extsize]} --nomodel -g {param[species]}  {param[treat_opt]} -n {param[description]}
                awk  \'{{OFS="\\t";if($5>1000){{$5=1000}}; print $0}}\' {output[peaks]} > {output[narrow]}

                fi
                eap_narrowPeak_to_bigBed {param[chrom]} {param[as]} {output[narrow]}  {output[bb]}
                wc -l {output[narrow]} > {output[qc]}
                """,
                tool="ls",
                input={"treat":target + "_final.bam"},
                output={"peaks": target + "_peaks.narrowPeak",
                        "narrow": target + ".narrowPeak",
                        "qc": target + ".narrowPeak.qc",
                        "bb": target + ".narrowPeak.bigBed",
                        "summit": target + "_summits.bed",
                        "treat_bdg": target + "_treat_pileup.bdg",
                        "peaks_xls": target + "_peaks.xls",
                        "control_bdg": target + "_control_lambda.bdg"},
                param={"description": target,
                       "tool": conf.get("tool", "peak_calling"),
                       "chrom": conf.get(conf.species, "chrom_len"),
                       "as": narrow,
                       "keep_dup": "all", "extsize": 50, "species": "hs" if conf.species in ["hg19","hg38"] else "mm",
                       "fdr": 0.01},
                name="macs2_callpeak_rep"))

        macs2_on_rep.param["treat_opt"] = " -t " + macs2_on_rep.input["treat"]
        macs2_on_rep.update(param=conf.items("macs2"))
        if conf.pe:
            macs2_on_rep.param["format"] = " -f BAMPE  "   ## pair end mode, not compatible with SAM files built-in sampling
        else:
            macs2_on_rep.param["format"] = " -f BAM "
        spot = attach_back(workflow,
            ShellCommand("""
            edwBamStats -sampleBamSize=5000000 -u4mSize=5000000 -sampleBam={input[tags]}.5M.bam {input[tags]} tmp.ra
            {tool} {input[tags]}.5M.bam {input[peaks]} {output}
            """,
            tool = "dac_macs2_spot.sh",
            input = {"tags": target + "_final_nochrm.bam",
                     "peaks": target + "_peaks.narrowPeak"},
            output = target + "_spot_nochrm_5M.qc",
            name = "macs2 spot"))
        bdg2bw = attach_back(workflow,
                             ShellCommand(
                                 "{tool} {input} {param[chromsize]} {output}",
                                 tool = "bdg2bw",
                                 input = target + "_treat_pileup.bdg",
                                 output = target + ".bigWig",
                                 param = {"chromsize": conf.get(conf.species, "chrom_len")}))
    have_treat_reps = len(conf.treatment_pairs) >= 2 ## replicates
    if have_treat_reps:
        eval_reps(workflow, conf, tex)
        catsam = attach_back(workflow, ShellCommand(
            "{tool} cat {param[bams]} > {output[bam]}",
            tool = "samtools",
            input ={"bams": [ target + "_final.bam" for target in conf.treatment_targets]},
            output = {"bam": conf.prefix + "_pool.bam"}))
        catsam.param.update(bams=' '.join(catsam.input["bams"]))
        if conf.species in ["hg19","hg38"]:
            species = "hs"
        else :
            species = "mm"

        macs_merge = macs2_on_rep.clone
        macs_merge.name = " macs2 merge "
        macs_merge.param = {"treat_opt": " -t " + conf.prefix + "_pool.bam",
                            "chrom": conf.get(conf.species, "chrom_len"),
                            "as": narrow,
                            "tool": conf.get("tool", "peak_calling"),
                            "description": conf.prefix, "keep_dup": "all", "extsize": 50, "species": species, "fdr": 0.01}
        if conf.pe:
            macs_merge.param["format"] = " -f BAMPE  "   ## pair end mode, not compatible with SAM files built-in sampling
        else:
            macs_merge.param["format"] = " -f BAM "

        macs_merge.output={"peaks": conf.prefix + "_peaks.narrowPeak",
                           "narrow": conf.prefix + ".narrowPeak",
                           "bb": conf.prefix + "narrowPeak.bigBed",
                           "summit": conf.prefix + "_summits.bed",
                           "treat_bdg": conf.prefix + "_treat_pileup.bdg",
                           "peaks_xls": conf.prefix + "_peaks.xls",
                           "qc": conf.prefix + ".narrowPeak.qc",
                           "control_bdg": conf.prefix + "_control_lambda.bdg"}
        attach_back(workflow, macs_merge)
        bdg_merge = bdg2bw.clone
        bdg_merge.input = conf.prefix + "_treat_pileup.bdg"
        bdg_merge.output = conf.prefix + ".bigWig"
        attach_back(workflow, bdg_merge)

