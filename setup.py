from distutils.core import setup

setup(
    name='DNP',
    version='0.1',
    packages=['dnp'],
    url='',
    license='MIT',
    author='Qian Qin',
    author_email='qinqianhappy@gmail.com',
    description='DNase I seq QC pipeline',
    scripts = ["dnp/DNP.py", "dnp/pipeline-scripts/tags.sh", "dnp/pipeline-scripts/script-tokenizer.py", "dnp/pipeline-scripts/runhotspot"],
    requires=["samflow", "jinja2", "argparse"],
    package_data = {"dnp" : ["static/*", "pipeline-scripts/*"]}
)
