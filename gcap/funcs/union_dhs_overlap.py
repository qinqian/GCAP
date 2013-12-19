#########################################################
#
# union DHS evaluation
# add cutting bias borrowed from Tarela
#
#########################################################
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.helpers import *
from pkg_resources import resource_filename

def union_DHS_overlap(workflow, conf, tex):
    """ use hotspot d narrow peaks for DHS evaluation """
    for i, target in enumerate(conf.treatment_targets):
        if conf.peakcalltool == "hotspot":
            peaks = conf.hotspot_reps_final_5M_prefix[i] + ".fdr0.01.pks.bed"
            output = target + "_DHS_overlap_peaks_bed"
        elif conf.peakcalltool == "macs2":
            peaks = target + "_5M_macs2_velcro_non_overlap_peaks.bed.final"
            output = target + "_macs2_DHS_overlap_peaks_bed"

        attach_back(workflow,
            ShellCommand(
                "{tool} -wa -u  \
                -a {input[pks_spot_bed]} -b {input[DHS_peaks_bed]} > {output}",
                tool="intersectBed",
                input={"pks_spot_bed": peaks,
                       "DHS_peaks_bed": conf.get("lib", "dhs")},
                output = output,
                name = "Write out DHS overlap BED"))

        if conf.get("cut_bias", "run").strip().upper() == "T":
            ## overall mode with make_bg_bwDHS_twobit.py
            if conf.get("cut_bias", "mode") == "overall":
                attach_back(workflow,
                            ShellCommand(
                                """
                                {tool} bamtobed -i {input[tag]} > {output[tags]}
                                awk '{{if ($6=="+") print $0}}' {output[tags]} > {output[plus_tag]}
                                awk '{{if ($6=="-") print $0}}' {output[tags]} > {output[minus_tag]}
                                macs2 pileup -i {output[plus_tag]} -o {output[plus_bdg]} -f BED --extsize 1
                                macs2 pileup -i {output[minus_tag]} -o {output[minus_bdg]} -f BED --extsize 1
                                bedClip {output[plus_bdg]} {param[chr]} {output[plus_bdg]}.clip
                                bedGraphToBigWig -unc {output[plus_bdg]}.clip {param[chr]} {output[bigwig_plus]}
                                bedClip {output[minus_bdg]} {param[chr]} {output[minus_bdg]}.clip
                                bedGraphToBigWig -unc {output[minus_bdg]}.clip {param[chr]} {output[bigwig_minus]}
                                rm -f {output[plus_bdg]}.clip
                                rm -f {output[minus_bdg]}.clip
                                """,
                                tool = "bedtools",
                                input = {"tag": target + ".bam"},
                                output = {"tags": target + "_cut_input.bed", "bigwig_plus": target + "_plus.bw",
                                          "bigwig_minus": target + "_minus.bw",
                                          "plus_tag": target + "cut_plus.bed", "minus_tag": target + "cut_minus.bed",
                                          "plus_bdg": target + "plus.bdg", "minus_bdg": target + "minus.bdg"},
                                name = " extract 1bp pileup from bam ",
                                param = {"chr": conf.get("lib", "chrom_len")}))


                attach_back(workflow,
                            ShellCommand(
                                """
                                ## use whole genome for atac seq
                                {tool} {param[script]} -p {input[bigwig_plus]} -n {input[bigwig_minus]} --genome={param[genome]} --sc {param[nmer_cut]} \
                                --left={param[mer]} --right={param[mer]} -o {output[matrix]} -i {input[peaks]}""",
                                tool = "python2.7",
                                input = {"bigwig_minus": target + "_minus.bw", "bigwig_plus": target + "_plus.bw",
                                         "peaks": conf.get("lib", "dhs")},
                                output = {"matrix": target + "_cut_bias.xls"},
                                param = {"script": resource_filename("gcap", "pipeline-scripts/make_bg_bwDHS_twobit.py"),
                                         "genome": conf.get("cut_bias", "seq_2bit"), ## "peaks": conf.get("lib", "dhs"),
                                         "nmer_cut": 10,
                                         "mer": 3}))

                ## Shawn's version of footprint
                ### 1. overlap with your DHS site , bedtools required
                attach_back(workflow,
                        ShellCommand(
                            """
                            ## input is chip seq peaks, openchr is dnase, atac or other open chromatin techniques peaks
                            ## support for macs2 peaks calling mode only
                            cut -f 1,2,3,4,9 {input[dhs]} > {input[dhs]}.bed
                            intersectBed -a {input[motif]} -b {input[dhs]}.bed -wa -u > {output[motif_dhs]} """,
                            tool = "intersectBed",
                            input = {"motif": conf.get("footprint", "motif"),
                                     "dhs": target + "_macs2_all_peaks.narrowPeak"},
                            output = {"motif_dhs": target + "_motif_dhs.bed"}))

                ### 2. scan for plotting data use python script, Scan_6c_matrix_cv_chip.py
                attach_back(workflow,
                        ShellCommand(
                            """
                            {tool} {param[script]} -i {input[motif_bed]} -o {output[matrix]} -b {input[cut_bias]} -p {input[bigwig_plus]} -n {input[bigwig_minus]} -c {input[conservation]} --chipbw {input[chip]}
                            """,
                            tool = "python2.7",
                            input = {"bigwig_plus": target + "_plus.bw",
                                     "bigwig_minus": target + "_minus.bw",
                                     "cut_bias": target + "_cut_bias.xls",
                                     "chip": conf.get("footprint", "chip"),
                                     "motif_bed": target + "_motif_dhs.bed", "conservation": conf.get("footprint", "conservation")},
                            output = {"matrix": target + "_footprint.xls"},
                            param = {"script": "Scan_6c_matrix_cv_chip.py"}))

                attach_back(workflow,
                    ShellCommand(
                        """
                        Rscript plot6c_cv_chip_10k.r {input[matrix]} {output}
                        """,
                        tool = "Rscript",
                        input = {"matrix": target + "_footprint.xls"},
                        output = target + "_footprint")) 
                ## pyDNase involved in footprint detection
                ## pyDNase
                ## import pyDNase

            ## strand mode for cutting bias
            ## cutting bias in stranded mode, a bit strange for ATAC-seq
            if conf.get("cut_bias", "mode") == "strand":
                if "bed" in conf.seq_type:
                    cut = attach_back(workflow,
                        ShellCommand(
                            "{tool} {param[script]} -p {input[peaks]} -s {input[bit]} -o {output[matrix]} -f {param[mer]} -t {input[tag]}",
                            tool = "python2.7",
                            input = {"peaks": conf.get("lib", "dhs"), "bit": conf.get("cut_bias", "seq_2bit"),
                                     "tag": target + "_all.bed"},
                            output = {"matrix": conf.treatment_targets[i] + "_cut_bias.xls"},
                            param = {"mer": 3, "script": resource_filename("gcap", "pipeline-scripts/twobit_seqbias.py")},
                            name = "cutting bias"))
                    cut.update(param = conf.items("cut_bias"))
                else:
                    cut = attach_back(workflow,
                        ShellCommand(
                            "bamToBed -i {input[tag]} > {output[tag]} && {tool} {param[script]} -p {input[peaks]} -s {input[bit]} -o {output[matrix]} -f {param[mer]} -t {output[tag]}",
                            tool = "python2.7",
                            input = {"peaks": conf.get("lib", "dhs"), "bit": conf.get("cut_bias", "seq_2bit"),
                                     "tag": target + ".bam"},
                            output = {"tag": target + "_cut_input.bed", "matrix": conf.treatment_targets[i] + "_cut_bias.xls"},
                            param = {"mer": 3, "script": resource_filename("gcap", "pipeline-scripts/twobit_seqbias.py")},
                            name = "cutting bias"))
                    cut.update(param = conf.items("cut_bias"))


    attach_back(workflow, PythonCommand(
        stat_dhs,
        input={"dhs_peaks": [ target + "_DHS_overlap_peaks_bed" if conf.peakcalltool == "hotspot" else target + "_macs2_DHS_overlap_peaks_bed"
                              for target in conf.treatment_targets ],
               "pks_spot_bed": [ conf.hotspot_reps_final_5M_prefix[i] + ".fdr0.01.pks.bed" if conf.peakcalltool == "hotspot" else target + "_5M_macs2_velcro_non_overlap_peaks.bed.final"
                                 for i, target in enumerate(conf.treatment_targets) ]},
        output = {"json": conf.json_prefix + "_dhs.json"},
        param= {"samples":conf.treatment_bases},
        name="DHS summary"))
    attach_back(workflow, PythonCommand(
        DHS_doc,
        input = {"tex": tex, "json": conf.json_prefix + "_dhs.json"},
        output = {"latex": conf.latex_prefix + "_dhs.tex"},
        param = {"reps": len(conf.treatment_pairs), "samples": conf.treatment_bases}))


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

def DHS_doc(input = {"json": "", "tex": ""}, output = {"latex": ""}, param = {"reps": "", "samples": ""}):
    data = json_load(input["json"])["stat"]
    dhs = []

    for s in param["samples"]:
        dhs.append(decimal_to_latex_percent(data[s]["dhspercentage"]))

    DHS_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "DHS",
        param = {"section_name": "DHS",
                 "render_dump": output["latex"],
                 "DHS": dhs,
                 "reps": param["reps"]})
    template_dump(DHS_latex)
