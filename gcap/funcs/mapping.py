#########################################################
#
# Read mapping QC
#
#########################################################

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back
from pkg_resources import resource_filename
from gcap.funcs.helpers import *

## default bwa command line for pair end and single end DNase
def _bwa(workflow, conf):
    """
    incorpate ENCODE ChIP-seq alignment parameters
    """
    for raw, target in conf.treatment_pairs:
        param = {"threads": conf.threads,
                 "index":conf.get(conf.species, "genome_index"),
                 "prefix": target + "_raw_sorted",
                 "qc2": target + "_rawbam_stats.qc"}

        if conf.pe:
            bwa = attach_back(workflow, ShellCommand(
                "{tool} {param[threads]} {param[index]} {input[fastq][0]} {input[fastq][1]} {output[bam]} {output[qc]} {param[prefix]} {param[qc2]}",
                tool = "eap_run_bwa_pe",
                input = {"fastq": raw},
                output = {"bam": target + "_raw_sorted.bam", "qc": target + "_rawbam.qc"},
                param = param,
                name = "pair end mapping"))
        else:
            bwa = attach_back(workflow, ShellCommand(
                "{tool} {param[threads]} {param[index]} {input[fastq]} {output[bam]} {output[qc]} {param[prefix]} {param[qc2]}",
                tool = "eap_run_bwa_se",
                input = {"fastq": raw},
                output = {"bam": target + "_raw_sorted.bam", "qc": target + "_rawbam.qc"},
                param = param,
                name = "single end mapping"))
        bwa.update(param = conf.items("bwa"))


## reads mapping by customized mapping tool, calculate autosome(exclude M, X, Y) mapping ratio as we discussed before
def read_mapping(workflow, conf, tex):
    """
    reads mapping for all reads
    and convert sam to bam, sort
    Sam or Bam input, only do statistics
    """
    _bwa(workflow, conf)

    # attach_back(workflow, PythonCommand(
    #     reads_doc,
    #     input = {"tex": tex, 
    #              "json_autosome": conf.json_prefix + "_um_autosome.json"},
    #     output = {"raw": conf.latex_prefix + "raw.tex",
    #               "mapping": conf.latex_prefix + "_mapping.tex"},
    #     param = {"seq_type": conf.pe,
    #              "reps": len(conf.treatment_pairs),
    #              "samples": conf.treatment_bases}))


        ## convert to BAM file
        # for n, target in enumerate(conf.treatment_targets):
        #     if not conf.seq_type.startswith("bam"):
        #         attach_back(workflow,
        #             ShellCommand(
        #                 "{tool} view -bq {param[map_quality]} -t {input[chrom_len]} {input[sam]} -o {output[bam]}",
        #                 tool="samtools",
        #                 input={"sam": target + "_all.sam", "chrom_len": conf.get_path("lib", "chrom_len")},
        #                 output={"bam": target + ".bam"},
        #                 param = {"map_quality": 1} ## to filter unreliable mapping
        #             ))
        #         workflow.update(param=conf.items("lib"))
        #     else:
        #         ## filtered BAM input Files do not need to convert
        #         link = attach_back(workflow, ShellCommand(
        #             ## filter or not, it's up to the user
        #             "{tool} -sf {input[bam]} {output[bam]}",
        #             tool = "ln",
        #             input = {"bam": conf.treatment_bam[n]},
        #             output = {"bam": target + ".bam"},   ## actually, may not be sorted bam files
        #             param = None))
