__author__ = 'qinqianhappy'

#########################################################
#
# Integrate all latex report to pdf
#
#########################################################
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from gcap.funcs.helpers import *

def report(workflow, conf, tex):
    collector = []

    ## begin doc and table header
    attach_back(workflow, PythonCommand(
        begin_doc,
        input = tex,
        output = {"begin": conf.latex_prefix + "_begin.tex",
                  "header": conf.latex_prefix + "_header.tex"},
        param = {"rep": len(conf.treatment_pairs), "id": underline_to_space(conf.id)}))
    # end doc
    attach_back(workflow, PythonCommand(
        end_doc,
        input = {"tex": tex},
        output = {"table_end": conf.latex_prefix + "_table_end.tex",
                  "doc_end": conf.latex_prefix + "_doc_end.tex"}))

    collector.append(conf.latex_prefix + "_begin.tex")
    collector.append(conf.latex_prefix + "_header.tex")
    if not conf.seq_type.startswith("bed"):
        collector.append(conf.latex_prefix + "raw.tex")
        collector.append(conf.latex_prefix + "len.tex")
        collector.append(conf.latex_prefix + "seq_quality.tex")
        collector.append(conf.latex_prefix + "_mapping.tex")
        collector.append(conf.latex_prefix + "_nsc.tex")
        collector.append(conf.latex_prefix + "_rsc.tex")

    collector.append(conf.latex_prefix + "_frag.tex")
    collector.append(conf.latex_prefix + "redun.tex")
    collector.append(conf.latex_prefix + "_peaks.tex")
    collector.append(conf.latex_prefix + "_spot.tex")
    if len(conf.treatment_pairs) >= 2:
        collector.append(conf.latex_prefix + "_reps.tex")
        ## not needed for ENCODE
    #    collector.append(conf.latex_prefix + "_promotor.tex")
    #    collector.append(conf.latex_prefix + "_conserv.tex")
    collector.append(conf.latex_prefix + "_dhs.tex")
    collector.append(conf.latex_prefix + "_table_end.tex")
    collector.append(conf.latex_prefix + "_doc_end.tex")

    integrate_doc = attach_back(workflow, ShellCommand(
        "cat {param[input]} > {output[tex]} \
         && {tool} -output-directory {output[dir]} -jobname={param[name]} {output[tex]}\
         && {tool} -output-directory {output[dir]} -jobname={param[name]} {output[tex]}",
        tool="pdflatex",
        input= collector[0],
        # output[pdf] should use "conf.prefix" to have the absolute path
        output={"tex": conf.prefix + "_report.tex", "dir": conf.target_dir, "pdf": conf.prefix + "_report.pdf"},
        # param[name] should use "conf.id" to avoid using absolute path
        param={"name": conf.id + "_report"},
        name="report"))
    integrate_doc.param.update({"input": " ".join(collector)})
