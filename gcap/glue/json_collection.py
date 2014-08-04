from gcap.funcs.helpers import *


def json_collect(conf):
    for target in conf.treatment_targets:
        qc = conf.target + ""
            ## not skip fastqc and do exist json file
        if exist(js):
            json_load("")
    json_dump(output["json"])
    return

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
# def _peaks_calling_latex(workflow, conf, tex):
#     ## peaks number json
#     if conf.peakcalltool == "hotspot":
#         if len(conf.treatment_pairs) >= 2:
#             ## use hotspot b to evaluate 5M reads peaks number, replicates consistency
#             ## use peaks d as narrow peaks to evaluate all reads peaks number
#             attach_back(workflow, PythonCommand(
#                 stat_peaks,
#                 input = {"peaks": {"5M_spot": [ conf.hotspot_reps_final_5M_prefix[i] + ".hot.bed" for i in range(len(conf.treatment_pairs)) ],
#                                    "all_peaks": [ conf.hotspot_reps_final_prefix[i] + ".fdr0.01.pks.bed" for i in range(len(conf.treatment_pairs)) ],
#                                    "combo": conf.hotspot_merge_final_prefix + "_merge_all.fdr0.01.pks.bed" }},
#                 output = {"json": conf.json_prefix + "_peaks.json"},
#                 param = {"tool": conf.peakcalltool, "samples": conf.treatment_bases},
#                 name = "json peaks and hotspot"))
#         else:
#             attach_back(workflow, PythonCommand(
#                 stat_peaks,
#                 input = {"peaks": {"5M_spot": [ conf.hotspot_reps_final_5M_prefix[i] + ".hot.bed" for i in range(len(conf.treatment_pairs)) ],
#                                    "all_peaks": [ conf.hotspot_reps_final_prefix[i] + ".fdr0.01.pks.bed" for i in range(len(conf.treatment_pairs)) ]}},
#                 output = {"json": conf.json_prefix + "_peaks.json"},
#                 param = {"tool": conf.peakcalltool, "samples": conf.treatment_bases},
#                 name = "json peaks and hotspot"))

#     elif conf.peakcalltool == "macs2":
#         if len(conf.treatment_pairs) >= 2:
#             attach_back(workflow, PythonCommand(
#                 stat_peaks,
#                 input = {"peaks": {"all_peaks": [ target + "_macs2_all_peaks.encodePeak" for target in conf.treatment_targets ],
#                                    "5M_spot": [ target + "_5M_macs2_peaks.encodePeak" for target in conf.treatment_targets ],
#                                    "combo": conf.prefix + "_peaks.encodePeak"}},
#                 output = {"json": conf.json_prefix + "_peaks.json"},
#                 param = {"tool": conf.peakcalltool, "samples": conf.treatment_bases}))
#         else:
#             attach_back(workflow, PythonCommand(
#                 stat_peaks,
#                 input = {"peaks": {"all_peaks": [ target + "_macs2_all_peaks.encodePeak" for target in conf.treatment_targets ],
#                                    "5M_spot": [ target + "_5M_macs2_peaks.encodePeak" for target in conf.treatment_targets ]}},
#                 output = {"json": conf.json_prefix + "_peaks.json"},
#                 param = {"tool": conf.peakcalltool, "samples": conf.treatment_bases}))

#     ## peaks number latex
#     attach_back(workflow, PythonCommand(
#         peaks_doc,
#         input = {"tex": tex, "json": conf.json_prefix + "_peaks.json"},
#         output = {"latex": conf.latex_prefix + "_peaks.tex"},
#         param = {"reps": len(conf.treatment_pairs), "samples": conf.treatment_bases},
#         name = "peaks number json"))

#     if conf.peakcalltool == "hotspot":
#         spot = [ target + "_5M_sort.spot.out" for target in conf.treatment_targets ]
#     elif conf.peakcalltool == "macs2":
#         spot = [ target + "_5M_macs2_peaks.encodePeak" + ".spot.out" for target in conf.treatment_targets ]

#     ## 5M reads, calculate the merged two passes and peaks regions number
#     ## 5M reads for macs2 optionally
#     attach_back(workflow, PythonCommand(
#         stat_spot_on_replicates,
#         input = {"spot_files": spot},
#         output = {"json": conf.json_prefix + "_sample_spot_5M.json"},
#         name = "json spot",
#         param = {"samples": conf.treatment_bases}))

#     attach_back(workflow, PythonCommand(
#         spot_doc,
#         input = {"json": conf.json_prefix + "_sample_spot_5M.json",
#                  "tex": tex},
#         output = {"latex": conf.latex_prefix + "_spot.tex"},
#         name = "doc spot 5M",
#         param = {"samples": conf.treatment_bases, "reps": len(conf.treatment_pairs)}))
