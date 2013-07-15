from distutils.core import setup

setup(
    name='GCAP',
    version='0.1',
    packages=['gcap', "gcap.funcs"],
    url='',
    license='MIT',
    author='Qian Qin',
    author_email='qinqianhappy@gmail.com',
    description='DNase I seq QC pipeline',
    scripts = ["gcap/GCAP", "gcap/pipeline-scripts/bed_duplicates.sh", "gcap/pipeline-scripts/script-tokenizer.py",
               "gcap/pipeline-scripts/runhotspot", "gcap/pipeline-scripts/macs2_spot.sh"],
    requires=["samflow", "jinja2", "argparse"],
    package_data = {"gcap" : ["static/*", "pipeline-scripts/*"]}
)