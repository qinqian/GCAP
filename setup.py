from distutils.core import setup

setup(
    name='GCAP',
    version='0.1',
    packages=['gcap'],
    url='',
    license='MIT',
    author='Qian Qin',
    author_email='qinqianhappy@gmail.com',
    description='DNase I seq QC pipeline',
    scripts = ["gcap/GCAP", "gcap/pipeline-scripts/tags.sh", "gcap/pipeline-scripts/script-tokenizer.py", "gcap/pipeline-scripts/runhotspot", "gcap/pipeline-scripts/sampling_sam_by_num.sh"],
    requires=["samflow", "jinja2", "argparse"],
    package_data = {"gcap" : ["static/*", "pipeline-scripts/*"]}
)