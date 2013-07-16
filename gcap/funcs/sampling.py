#########################################################
#
# 1. C version sampling for fastq, BAM
# 2. python version sampling fastq, SAM, BED
#
#########################################################

#########################################################
# using fastqStatsAndSubsample  sampleBam
#########################################################
import random
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back

def sample_reads(workflow, conf, N, format):
    """
    get random N reads from Fastq, SAM, BAM
    """
    if format == "fastq":
        if conf.seq_type == "se":
            for n, target in enumerate(conf.treatment_targets):
                attach_back(workflow,
                    PythonCommand(single_end_fastq_sampling,
                        input = {"fastq": conf.treatment_raws[n]},
                        param = {"random_number": N},
                        output = {"fastq_sample": target + "_100k.fastq"}))
        elif conf.seq_type == "pe":
            for raw, target in conf.treatment_pairs_pe:
                attach_back(workflow,
                    PythonCommand(pair_end_fastq_sampling,
                    input={"fastq": raw},
                    output={"fastq_sample": [ i + "_100k.fastq" for i in target ]},
                    param={"random_number": N}))     ## use 100kb random reads
    elif format == "sam":
        suffix = "100k" if N == 100000 else "5M"
        for target, s in zip(conf.treatment_targets, conf.treatment_bases):
            ## built-in sampling for raw reads
            if conf.seq_type.startswith("bam") or conf.seq_type.startswith("sam"):
                attach_back(workflow, PythonCommand(
                    sampling,
                    input = {"sam": target + "_all.sam"},
                    output = {"sam_sample": target + "_%s.sam" % suffix},
                    param = {"random_number": N,
                             "se_or_pe": conf.seq_type.strip().split(",")[1].strip().lower()}))
            else:
                attach_back(workflow, PythonCommand(
                    sampling,
                    input = {"sam": target + "_all.sam"},
                    output = {"sam_sample": target + "_%s.sam" % suffix},
                    param = {"random_number": N,
                             "se_or_pe": conf.seq_type}))
    elif format == "bed":
        for target in conf.treatment_targets:
            ## macs2 randsample sampling mappable reads from BED reads files
            attach_back(workflow, ShellCommand(
                "{tool} randsample -t {input[bed]} -n {param[rand_number]} -o {output[bed]}",
                tool = "macs2",
                input = {"bed": target + "_all.bed"},
                output = {"bed": target + "_5M.bed"}, ## output bed file, to be uniform with other sampling methods
                param = {"rand_number": N}))

#########################################################
# python version sampling fastq, SAM, BED
#########################################################

def single_end_fastq_sampling(input = {"fastq": ""}, output = {"fastq_sample": ""}, param = {"random_number": ""}):
    num_lines = sum(1 for _ in open(input["fastq"])) / 4

    rand_nums = sorted([random.randint(0, num_lines - 1) for _ in range(param["random_number"])])

    fastq = open(input["fastq"],"rU")
    fastq_sample = open(output["fastq_sample"], "w")

    cur_num = - 1
    written = 0

    for rand_num in rand_nums:
        while cur_num < rand_num:
            for i in range(4):
                fastq.readline()
            cur_num += 1

        for i in range(4):
            fastq_sample.write(fastq.readline())
        cur_num += 1
        written += 1
    assert written == param["random_number"]

    fastq_sample.close()
    fastq.close()

def pair_end_fastq_sampling(input = {"fastq": ""}, output = {"fastq_sample": ""}, param = {"random_number": ""}):
    def write_random_records(fqa, fqb, suba, subb, N=100000):
        records = sum(1 for _ in open(fqa)) / 4
        rand_records = sorted([random.randint(0, records - 1) for _ in range(N)])

        fha, fhb = open(fqa),  open(fqb)
        suba, subb = open(suba, "w"), open(subb, "w")
        rec_no = - 1
        written = 0

        for rr in rand_records:
            while rec_no < rr:
                rec_no += 1
                for i in range(4): fha.readline()
                for i in range(4): fhb.readline()
            for i in range(4):
                suba.write(fha.readline())
                subb.write(fhb.readline())
            rec_no += 1
            written += 1
        assert written == N
        suba.close()
        subb.close()
        fha.close()
        fhb.close()
    write_random_records(input["fastq"][0], input["fastq"][1], output["fastq_sample"][0], output["fastq_sample"][1], param["random_number"])

def sampling(input = {"sam": ""}, output = {"sam_sample": ""}, param = {"random_number": "", "se_or_pe": ""}):
    """
    get 5M reads from SAM, BED
    """
    num_lines = sum(1 for _ in open(input["sam"]))
    header_num = 0
    header = []
    with open(input["sam"]) as f:
        for line in f:
            if line.startswith("@"):
                header_num += 1
                header.append(line)
            else:
                break
    print(header_num)
    if param["se_or_pe"] == "se":
        rand_nums = sorted([random.randint(header_num, num_lines - 1) for _ in range(param["random_number"])])
        print(len(rand_nums))
        cur_num = -1
        written = 0
        with open(output["sam_sample"], "w") as fout:
            with open(input["sam"], "rU") as fin:
                for rand_num in rand_nums:
                    while cur_num < rand_num:
                        cur_num+=1
                        data = fin.readline()
                        if data.startswith("@"):
                            fout.write(data)

                    fout.write(fin.readline())
                    cur_num += 1
                    written += 1
        assert  written == param["random_number"]
    elif param["se_or_pe"] == "pe":
        random_pe = []
        cur_num = -1
        written = 0
        ## large memory version 1, faster, with paired reads sampling
#            for _ in range(int(param["random_number"]/2)):
#                num = random.choice(range(header_num, num_lines, 2))
#                random_pe.append(num)
#                random_pe.append(num+1)
#            rand_nums = random_pe
#            data = open(input["sam"]).readlines()
#            with open(output["sam_sample"], 'w') as f:
#               f.write("".join(header))
#               for i in rand_nums:
#                   f.write(data[i])
        ## small memory version 2, slower and simpler
        ## TODO: suggested by XiuXiu, open file twice, sampling odds and even separately
        ## this would lead to paired reads
        data_range = range(header_num, num_lines, 2)
        for _ in range(int(param["random_number"]/2)):
            num = random.choice(data_range)
            random_pe.append(num)
            random_pe.append(num+1)

        rand_nums = sorted(random_pe)
        with open(output["sam_sample"], "w") as fout:
            fout.write("".join(header))
            with open(input["sam"], "rU") as fin:
                for rand_num in rand_nums:
                    while cur_num < rand_num:
                        cur_num+=1
                        fin.readline()
                    fout.write(fin.readline())
                    cur_num += 1
                    written += 1
            assert  written == param["random_number"]