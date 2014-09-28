from samflow.command import ShellCommand
from samflow.workflow import attach_back
import os

def sample_bam_stat(workflow, conf, tex):
    """ sample non chrm bam to 15M for NSC and PBC
    sample non chrm bam to 5M for spot
    """
    for i, target in enumerate(conf.treatment_targets):
        ## for PE, use name sorted in order to calculate PBC
        input_bam = target + "_name_sorted.bam" if conf.pe else target + "_final_nochrm.bam"
        attach_back(workflow, ShellCommand(
            "{tool} {input[namesorted]} {param[run_spp]} {output[bamstat]} {output[sppstat]}  {param[pe]} {output[pbc]}",
            tool = "eap_dnase_stats",
            input = {"namesorted": input_bam},
            output = {"bamstat": target + "_bam_stat.qc",  ## 15M
                      "sppstat": target + "_spp.qc",
                      "pbc": target + "_final_nochrm_15M_pbc.qc"},
            param = {"pe": "pe" if conf.pe else "se",
                     "run_spp": conf.get("tool", "spp")}))

        if not "macs" in conf.get("tool", "peak_calling"):

            attach_back(workflow, ShellCommand(
                "{tool} {input[bamwithoutchrm]} {param[genome]} {param[readsize]} {output[spot]} {param[hotspot_dir]} {param[hotspot_output]} {param[hotspot_tmp]} {param[spot_tmp]}",
                tool = "dac_spot", ## 5M
                input = {"bamwithoutchrm": target + "_final_nochrm.bam"},
                output = {"spot": target + "_spot_nochrm_5M.qc"},

                param = {"genome": conf.species,
                         "spot_tmp": conf.hotspot_reps_tmp_prefix[i] + "_final_nochrm.bam.5000000.spot.out",
                         "readsize": conf.readsize,
                         "hotspot_dir": conf.get("tool", "peak_calling"),
                         "hotspot_output": target + "_hotspot",
                         "hotspot_tmp": target + "_hotspot_tmp"}))
