#########################################################
#
# 1. C version sampling for fastq, BAM
#    using fastqStatsAndSubsample  sampleBam
# 2. python version sampling fastq, SAM, BED
#
#########################################################

from samflow.command import ShellCommand
from samflow.workflow import attach_back


def filter_bam(workflow, conf, tex):
    """ filter bam file by samtools and sample by ucsc app
    """

    for target in conf.treatment_targets:
        input = {"raw": target + "_raw_sorted.bam"}
        if conf.pe:
            name = "pair"
            tool = "dac_bam_pe_post_filter"
            param = {"mapq": 3, "namesortedbamprefix": target + "_name_sorted",
                     "finalprefix": target + "_final",
                     "qc2": target + "_filter_bam_stats.qc"}
            output = {"finalbam": target + "_final.bam",
                      "namesortedbam": target + "_name_sorted.bam",
                      "bamwithoutchrm": target + "_final_nochrm.bam",
                      "qc": target + "_filter_bam.qc"}

            attach_back(workflow, ShellCommand(
                "{tool} {input[raw]} {param[namesortedbamprefix]} {output[namesortedbam]} {param[finalprefix]} {output[finalbam]} {param[mapq]} {output[bamwithoutchrm]} {output[qc]} {param[qc2]}",
                tool=tool,
                input=input,
                output=output,
                param=param,
                name="%s end filtering" % name))
        else:
            name = "single"
            tool = "dac_bam_se_post_filter"
            param = {"mapq": 3,
                     "finalprefix": target + "_final",
                     "qc2": target + "_filter_bam_stats.qc"}
            output = {"finalbam": target + "_final.bam",
                      "bamwithoutchrm": target + "_final_nochrm.bam",
                      "qc": target + "_filter_bam.qc"}


            attach_back(workflow, ShellCommand(
                "{tool} {input[raw]} {output[finalbam]} {param[mapq]} {output[qc]} {output[bamwithoutchrm]} {param[finalprefix]} {param[qc2]}",
                tool=tool,
                input=input,
                output=output,
                param=param,
                name="%s end filtering" % name))
