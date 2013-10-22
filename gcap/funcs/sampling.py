#########################################################
#
# 1. C version sampling for fastq, BAM
#    using fastqStatsAndSubsample  sampleBam
# 2. python version sampling fastq, SAM, BED
#
#########################################################

from samflow.command import ShellCommand
from samflow.workflow import attach_back

def sample_reads(workflow, conf, N, format):
    """
    get random N reads from Fastq, SAM, BAM
    """
    if format == "fastq":
        if conf.seq_type == "se" or conf.seq_type.strip().split(',')[0].strip().lower() in ['bam', 'sam']:
            for n, target in enumerate(conf.treatment_targets):
                if conf.seq_type == "se":
                    data = conf.treatment_raws[n]
                if conf.seq_type.strip().split(',')[0].strip().lower() in ['bam', 'sam']:
                    data = target + "_100k.fastq.tmp"
                attach_back(workflow,
                    ShellCommand(
                        "{tool} -sampleSize={param[random_number]} {input} {output[stat]} {output[fastq]}",
                        tool = "fastqStatsAndSubsample",
                        input = data,
                        output = {"fastq": target + "_100k.fastq",
                                  "stat": target + "_100k.seq"},
                        param = {"random_number": N}))

        elif conf.seq_type == "pe":
            for raw, target in conf.treatment_pairs_pe:
                attach_back(workflow,
                    ShellCommand(
                        "{tool} -sampleSize={param[random_number]} {input[fastq][0]} {output[stat][0]} {output[fastq][0]} && \
                        {tool} -sampleSize={param[random_number]} {input[fastq][1]} {output[stat][1]} {output[fastq][1]}",
                        tool = "fastqStatsAndSubsample",
                        input = {"fastq": raw},
                        output = {"fastq": [ i + "_100k.fastq" for i in target ],
                                  "stat": [ i + "_100k.seq" for i in target ]},
                        param = {"random_number": N}))

    elif format == "sam":
        suffix = "100k" if N == 100000 or N == 110000 else "5M"
        for target, s in zip(conf.treatment_targets, conf.treatment_bases):
            ## samtools sampling version, sampling by ratio, always more than 100k or 5M
            attach_back(workflow, ShellCommand(
                "a=$(wc -l {input[sam]} | cut -f 1 -d\" \") && b=$(echo \"scale=5;{param[random_number]}/$a\" | bc) \
                && echo $b && {tool} view -Shs $b {input[sam]} > {output[sam_sample_tmp]}",
                tool = "samtools",
                input = {"sam": target + "_all.sam"},
                output = {"sam_sample_tmp": target + "_%s.sam.tmp" % suffix},
                param = {"random_number": 1.05*N, ## get more than N number, in order to get top N reads
                         "se_or_pe": conf.seq_type}))
            ## get top 5M or 100k reads
            attach_back(workflow, ShellCommand(
                "a=$({tool} view -SH {input[sam_sample_tmp]} | wc -l) && a=$(($a+{param[random_number]})) && head -n $a {input[sam_sample_tmp]} > {output[sam_sample]}",
                tool = "samtools",
                input = {"sam_sample_tmp": target + "_%s.sam.tmp" % suffix},
                output = {"sam_sample": target + "_%s.sam" % suffix},
                param = {"random_number": N}))

    elif format == "bam":
        suffix = "100k" if N == 100000 or N == 110000 else "5M"
        for target, s in zip(conf.treatment_targets, conf.treatment_bases):
            attach_back(workflow, ShellCommand(
                "a=$({tool} flagstat {input[bam]} | head -1 | cut -f 1 -d\" \") && b=$(echo \"scale=5;{param[random_number]}/$a\" | bc) \
                && echo $b && {tool} view -bhs  $b {input[bam]} > {output[bam_sample]}",
                tool = "samtools",
                input = {"bam": target + "_all.bam"},
                output = {"bam_sample": target + "_%s.bam.tmp" % suffix},
                param = {"random_number": 1.05*N,
                         "se_or_pe": conf.seq_type}))

            ## get top 5M or 100k reads
            attach_back(workflow, ShellCommand(
                "a=$({tool} view -H {input[bam_sample_tmp]} | wc -l) && a=$(($a+{param[random_number]})) && \
                samtools view -h {input[bam_sample_tmp]} | head -n $a > {output[sam_sample]}",
                tool = "samtools",
                input = {"bam_sample_tmp": target + "_%s.bam.tmp" % suffix},
                output = {"sam_sample": target + "_%s.sam" % suffix},
                param = {"random_number": N}))

    elif format == "bed":
        for i, target in enumerate(conf.treatment_targets):
            ## macs2 randsample sampling mappable reads from BED reads files
            attach_back(workflow, ShellCommand(
                "{tool} randsample -t {input[bed]} -n {param[rand_number]} -o {output[bed]}",
                tool = "macs2",
                input = {"bed": target + "_all.bed"},
                output = {"bed": target + "_5M.bed"}, ## output bed file, to be uniform with other sampling methods
                param = {"rand_number": N}))