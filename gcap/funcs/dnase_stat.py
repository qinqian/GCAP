from samflow.command import ShellCommand
from samflow.workflow import attach_back

def sample_bam_stat(workflow, conf, tex):
    """ sample non chrm bam to 15M for NSC and PBC
    sample non chrm bam to 5M for spot
    """
    for target in conf.treatment_targets:
        attach_back(workflow, ShellCommand(
            "{tool} {input[namesorted]} {param[run_spp]} {output[bamstat]} {output[sppstat]}  {param[pe]} {output[pbc]}",
            tool = "eap_dnase_stats",
            input = {"namesorted": target + "_name_sorted.bam"},
            output = {"bamstat": target + "_bam_stat.qc",  ## 15M
                      "sppstat": target + "_spp.qc",
                      "pbc": target + "_final_nochrm_15M_pbc.qc"},
            param = {"pe": "pe" if conf.pe else "se",
                     "run_spp": conf.get("tool", "spp")}))

        attach_back(workflow, ShellCommand(
            "{tool} {input[bamwithoutchrm]} {param[genome]} {param[readsize]} {output[spot]} {param[hotspot_dir]} {param[hotspot_output]} {param[hotspot_tmp]}",
            tool = "dac_spot", ## 5M
            input = {"bamwithoutchrm": target + "_final_nochrm.bam"},
            output = {"spot": target + "_spot_nochrm_5M.qc"},
            param = {"genome": conf.species,
                     "readsize": conf.readsize,
                     "hotspot_dir": conf.get("tool", "peak_calling"),
                     "hotspot_output": target + "_hotspot",
                     "hotspot_tmp": target + "hotspot_tmp"}))
