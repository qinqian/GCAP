from gcap.funcs.helpers import *


def json_collect(input="",output="",param=""):
    json_dump(output["json"])
    return

def stat_strand_cor(input = {"metric": ""}, output = {"json": ""}, param= {"samples": ""}):
    """RSC NSC: COL9 COL10, COL3: estFragLen"""
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}

    for i, j in zip(input["metric"], param["samples"]):
        data = open(i).readlines()
        l = data[0].strip().split("\t")
        json_dict["stat"][j] = {}
        json_dict["stat"][j]["RSC"] = l[8]
        json_dict["stat"][j]["NSC"] = l[9]
        json_dict["stat"][j]["Frag"] = l[2]
    json_dump(json_dict)

## parse fastqStatsAndSubsample result
def stat_fastqStat(input = {"seq": ""}, output = {"json": ""}, param = {"samples": "", "seq_type": ""}):
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    if param["seq_type"] == "pe":
        for i, s in zip(input['seq'], param['samples']):
            with open(i[0]) as f:
                data = f.readlines()
                json_dict["stat"][s + "_pair1"] = {}
                json_dict["stat"][s + "_pair1"]['count'] = data[0].strip().split()[1]
                json_dict["stat"][s + "_pair1"]['basecount'] = data[1].strip().split()[1]
                json_dict["stat"][s + "_pair1"]['quality'] = round(float(data[8].strip().split()[1]), 1)
                json_dict["stat"][s + "_pair1"]['std'] = data[9].strip().split()[1]
                json_dict["stat"][s + "_pair1"]['len'] = data[4].strip().split()[1]
            with open(i[1]) as f:
                data = f.readlines()
                json_dict["stat"][s + "_pair2"] = {}
                json_dict["stat"][s + "_pair2"]['count'] = data[0].strip().split()[1]
                json_dict["stat"][s + "_pair2"]['basecount'] = data[1].strip().split()[1]
                json_dict["stat"][s + "_pair2"]['quality'] = round(float(data[8].strip().split()[1]), 1)
                json_dict["stat"][s + "_pair2"]['std'] = data[9].strip().split()[1]
                json_dict["stat"][s + "_pair2"]['len'] = data[4].strip().split()[1]
    elif param["seq_type"] == "se":
        for i, s in zip(input['seq'], param['samples']):
            with open(i) as f:
                data = f.readlines()
                json_dict["stat"][s] = {}
                json_dict["stat"][s]['count'] = data[0].strip().split()[1]
                json_dict["stat"][s]['basecount'] = data[1].strip().split()[1]
                json_dict["stat"][s]['quality'] = data[8].strip().split()[1]
                json_dict["stat"][s]['std'] = data[9].strip().split()[1]
                json_dict["stat"][s]['len'] = data[4].strip().split()[1]
    json_dump(json_dict)

        # attach_back(workflow, PythonCommand(stat_fastqStat,
        #                                     input = {"seq": [ t + "_100k.seq" for t in conf.treatment_targets ]},
        #                                     output = {"json": conf.json_prefix + "_seq_quality.json"},
        #                                     param = {"samples": conf.treatment_bases,
        #                                              "seq_type": conf.pe}))
        # attach_back(workflow, PythonCommand(
        #     seq_quality_doc,
        #     input = {"tex": tex, "json": conf.json_prefix + "_seq_quality.json"},
        #     output = {"seq": conf.latex_prefix + "seq_quality.tex", "len": conf.latex_prefix + "len.tex"},
        #     param = {"seq_type": conf.seq_type, "reps": len(conf.treatment_pairs), "se_samples": conf.treatment_bases}))

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
