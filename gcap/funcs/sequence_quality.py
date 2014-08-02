#############################################################
# The first step for GCAP QC, add multiple format support
#    input: fastq, use fastqStatsAndSubsample
#    output: sequence quality and subsampled reads
#############################################################

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back

def read_quality(workflow, conf, tex):
    if conf.pe:
        for raw, target in conf.treatment_pairs_pe:
            attach_back(workflow,
                        ShellCommand(
                            "{tool} {input[fastq][0]} {input[fastq][1]} {output[stat][0]} {output[stat][1]}",
                            tool = "dac_pe_read_quality",
                            input = {"fastq": raw},
                            output = {"stat": [ i + "_read_quality.qc" for i in target ]}))

        # attach_back(workflow, PythonCommand(stat_fastqStat,
        #                                     input = {"seq": [ [ p + "_100k.seq" for p in target ] for target in conf.treatment_pair_data ]},
        #                                     output = {"json": conf.json_prefix + "_seq_quality.json"},
        #                                     param = {"samples": conf.treatment_bases, "seq_type": conf.pe}))
        # attach_back(workflow, PythonCommand(
        #     seq_quality_doc,
        #     input = {"tex": tex, "json": conf.json_prefix + "_seq_quality.json"},
        #     output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
        #     param = {"seq_type": conf.seq_type, "reps": len(conf.treatment_pairs),
        #              "pe_samples": conf.treatment_bases}))
    else:
        for raw, target in conf.treatment_pairs:
            sample_fq = {"stat": target + "_read_quality.qc"}
            attach_back(workflow,
                        ShellCommand(
                            "{tool} {input} {output[stat]}",
                            tool = "dac_se_read_quality",
                            input = raw,
                            output = sample_fq,
                            name = "100k read sequence quality and sequence length"))
