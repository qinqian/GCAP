#########################################################
#
# load jinja2 template and write out latex snippets
#   1. docs begining
#   2. docs ending
#
#########################################################

from configparser import ConfigParser
import json
import os
import re
import sys

from jinja2 import Environment, FileSystemLoader
from samflow.command import AbstractCommand
from samflow.workflow import attach_back

def json_dump(json_dict):
    json_file = json_dict["output"]["json"]
    with open(json_file, "w") as f:
        json.dump(json_dict, f, indent=4)
    return json_file

def json_load(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    return json_dict

env = Environment(loader = FileSystemLoader("/"),
    block_start_string = '\BLOCK{',
    block_end_string = '}',
    variable_start_string = '\VAR{',
    variable_end_string = '}',
    comment_start_string = '\#{',
    comment_end_string = '}',
    line_statement_prefix = '%-',
    line_comment_prefix = '%#',
    trim_blocks = True,
    autoescape = False,)

def count_in_million(x):
    if type(x) == int:
        return str(round(x/1000000, 1)) + "M"

def decimal_to_latex_percent(dec):
    if type(dec) == float:
        return str(round(dec*100, 1))    + "\%"

def underline_to_space(x):
    if type(x) == str:
        return x.replace("_", " ")
    return x

class JinjaTemplateCommand(AbstractCommand):
    def __init__(self, template, tool=None, param = {}, input=[], output=[], name = ""):
        AbstractCommand.__init__(self, template=template, tool=None, param = param, input=[], output=[], name = name)
        self.env = env
        self._t = self.env.get_template(self.template)

    def _execute(self):
        """ Real-run current command"""
        self.result = self._t.render(input = self.input, output = self.output, **self.param)
        return True

    def _simulate(self):
        """ Dry-run current command: Pretend to run but not invoke anything """
        print("Rendering Latex part %s" % self.name, self.template)
        return True

def template_dump(jinja_template):
    jinja_template.invoke()
    with open(jinja_template.param["render_dump"], "w") as f:
        f.write(jinja_template.result)


def begin_doc(input = "", output = {"begin": "", "header": ""}, param = {"rep": "", "id": ""}):
    ## begin document
    begin_latex = JinjaTemplateCommand(
        template = input,
        name = "begin",
        param = {"section_name": "begin",
                 "render_dump": output["begin"],
                 "dataset_id": param["id"]})
    template_dump(begin_latex)

    header_latex = JinjaTemplateCommand(
        template = input,
        name = "header",
        param = {"section_name": "header",
                 "render_dump": output["header"],
                 "reps": param["rep"]})
    template_dump(header_latex)

def end_doc(input = {"tex": ""}, output = {"table_end": "", "doc_end": ""}, param = {}):
    """
    table end and latex end
    """
    table_end_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "end",
        param = {"section_name": "table_end",
                 "render_dump": output["table_end"]})

    template_dump(table_end_latex)

    end_latex = JinjaTemplateCommand(
        template = input["tex"],
        name = "end",
        param = {"section_name": "end",
                 "render_dump": output["doc_end"]})
    template_dump(end_latex)
