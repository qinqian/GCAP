#!/usr/bin/env python
"""

## TODO: uniform json and better extraction, write into GCAP function

extract DNase data summary from GCAP pipeline
testid_conserv.json  testid_contam.json  testid_dhs.json  testid_frag.json  testid_peaks.json
testid_promotor.json  testid_redun.json  testid_rep1_fastqc.json  testid_rep2_fastqc.json  testid_reps.json  testid_sample_spot_5M.json  testid_um_autosome.json
"""

import os
import glob
import json

data = glob.glob("*/json")

def json_load(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    return json_dict

def json_dataset_parse(dataset):
    jsons = os.listdir(dataset)
    data = []
    for j in jsons:
        js = json_load(os.path.join(d, j))
        temp = []
        if "_um_autosome" in j:
            auto = []
            for i in js["stat"].values():
                auto.append(str(i.values()[2]))
            temp.append(";".join(auto))
            data+=temp

        if "_fastq" in j:
            fastq = []
            for i in js["stat"].values():
                fastq.append(str(i.values()[1]))
            temp.append(";".join(fastq))
            data+=temp

        if "_sample_spot_5M" in j:
            spot = []
            for i in js["stat"].values():
                spot.append(str(i.values()[1]))
            temp.append(";".join(spot))
            data+=temp

        if "redun" in j:
            temp.append(";".join(map(str,js["stat"].values())))
            data+=temp

        if "dhs" in j:
            dhs = []
            for i in js["stat"].values():
                dhs.append(str(i.values()[0]))
            temp.append(";".join(dhs))
            data+=temp

        if "_peaks" in j:
            temp.append(";".join(map(str, js["stat"]["all_peaks"].values())))
            temp.append(";".join(map(str,js["stat"]["5M_spot"].values())))
            data+=temp

        if "frag" in j:
            temp.append(";".join(map(str, js["stat"].values())))
            data+=temp

        if "NSC" in j:
            pass

        if "RSC" in j:
            pass

    return data

for d in data:
    if "real" in d or "encode" in d:
        continue
    else:
        result = json_dataset_parse(d)
        print(d, "\t".join(result))
