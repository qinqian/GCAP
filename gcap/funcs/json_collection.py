from gcap.funcs.helpers import *


def json_collect(input, output = {}, param = {}):
    conf = param
    collection = {"input": input, "output": output, "stat": {}}
    have_treat_reps = len(conf.treatment_pairs) >= 2 ## replicates

    for target,sample in zip(conf.treatment_targets, conf.treatment_bases):
        collection["stat"][sample] = {}

        if conf.pe:
            readquality1 = float(open(target + "pair1_read_quality.qc").readlines()[8].split()[1])
            readquality2 = float(open(target + "pair2_read_quality.qc").readlines()[8].split()[1])
            collection["stat"][sample]['ReadQuality'] = [readquality1, readquality2] ## pair 1,  pair 2
        else:
            readquality = float(open(target + "_read_quality.qc").readlines()[8].split()[1])
            collection["stat"][sample]['ReadQuality'] = [readquality] ## read

        raw = target + "_rawbam.qc"
        raw = open(raw).readlines()[0].split()[0]

        filter_bam = target + "_filter_bam.qc"
        filter_bam = open(filter_bam).readlines()[2].split()[0] ## mapping qc

        bamflag = target + "_bam_stat.qc" ## optional library complexity
        bamflag = open(bamflag).readlines()[13].strip().split()[1]

        bamstats = target + "_filter_bam_stats.qc" ## insert size
        bamstats = open(bamstats).readlines()[28].strip().split()[-1]

        narrowpeak = target + ".narrowPeak.qc"
        # broadpeak = target + ".broadPeak.qc"  ## peak number, total filtered uniquely mapped reads with redundancy

        phantom = target + "_spp.qc"
        phantom = open(phantom).read().strip().split()  ## NSC, 15M without chrM

        pbc = target + "_final_nochrm_15M_pbc.qc" ## no chrM, library complexity
        pbc = open(pbc).read().strip().split()

        collection["stat"][sample]["totalreads"] = int(raw)
        collection["stat"][sample]["mapped"] = int(filter_bam)

        collection["stat"][sample]["UniquePosRatio"] = float(bamflag)
        collection["stat"][sample]["PBC1"] = float(pbc[4])
        collection["stat"][sample]["FragmentSize"] = float(bamstats)
        # collection["stat"][sample]["RSC"] = float(phantom[8])
        collection["stat"][sample]["NSC"] = float(phantom[9])
        collection['stat'][sample]["narrowPeaknum"] = int(open(narrowpeak, 'rU').read().split()[0])
        # collection['stat'][sample]["broadPeaknum"] = int(open(broadpeak, 'rU').read().split()[0])
        spot = target + "_spot_nochrm_5M.qc"
        collection['stat'][sample]["spot"] = float(open(spot, 'rU').readlines()[1].split()[-1])
    #    for i, s in zip(input['seq'], param['samples']):
    #         with open(i[0]) as f:
    #             data = f.readlines()
    #             json_dict["stat"][s + "_pair1"] = {}
    #             json_dict["stat"][s + "_pair1"]['count'] = data[0].strip().split()[1]
    #             json_dict["stat"][s + "_pair1"]['basecount'] = data[1].strip().split()[1]
    #             json_dict["stat"][s + "_pair1"]['quality'] = round(float(data[8].strip().split()[1]), 1)
    #             json_dict["stat"][s + "_pair1"]['std'] = data[9].strip().split()[1]
    #             json_dict["stat"][s + "_pair1"]['len'] = data[4].strip().split()[1]
    #         with open(i[1]) as f:
    #             data = f.readlines()
    #             json_dict["stat"][s + "_pair2"] = {}
    #             json_dict["stat"][s + "_pair2"]['count'] = data[0].strip().split()[1]
    #             json_dict["stat"][s + "_pair2"]['basecount'] = data[1].strip().split()[1]
    #             json_dict["stat"][s + "_pair2"]['quality'] = round(float(data[8].strip().split()[1]), 1)
    #             json_dict["stat"][s + "_pair2"]['std'] = data[9].strip().split()[1]
    #             json_dict["stat"][s + "_pair2"]['len'] = data[4].strip().split()[1]
    collection['stat']['SeqType'] = "Pair End"if conf.pe else "Single End"
    if have_treat_reps:
        narrowpeak = conf.prefix + ".narrowPeak.qc"
        # broadpeak = conf.prefix + ".broadPeak.qc"
        cor = conf.prefix + "_cor.qc"
        overlap = conf.prefix + "_overlap.qc"
        collection['stat']['correlation_between_replicates'] = float(open(cor, 'rU').read().strip())
        collection['stat']['overlap_divided_by_union_bases'] = float(open(overlap, 'rU').readlines()[-1].strip().split()[1])
        collection['stat']['mergedNarrowPeak'] = int(open(narrowpeak, 'rU').read().split()[0])
        # collection['stat']['mergedBroadPeak'] = int(open(broadpeak, 'rU').read().split()[0])

    json_dump(collection)


def stat_dhs(input={"pks_spot_bed": "", "dhs_peaks": ""}, output={"json": ""},
             param={"samples": ""}):
    """
    overlap with union DHS
    """
    peaks_info = {}
    for b, d, s in zip(input["pks_spot_bed"], input["dhs_peaks"], param["samples"]):
        print(b, d, s)

        peaks_info[s] = {}
        peaks_info[s]["total"] = len(open(b, 'r').readlines())
        peaks_info[s]["dhs"] = len(open(d, 'r').readlines())
        peaks_info[s]['dhspercentage'] = peaks_info[s]["dhs"] / peaks_info[s]["total"]
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict["stat"] = peaks_info
    json_dump(json_dict)
