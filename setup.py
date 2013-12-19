"""
GCAP
-----

GCAP is a python module for wrapping different tools to quality control
DNase-seq and other open chromatin sequencing data.

GCAP is easy to start
`````````````````````

.. code:: python

    GCAP run -c config


Links
`````

* `website <https://github.com/qinqian/GCAP>`_
* `documentation <http://pythonhosted.org/GCAP/>`_


"""

from setuptools import setup

setup(
    name='GCAP',
    version='1.0',
    packages=['gcap', "gcap.funcs"],
    url='https://github.com/qinqian/GCAP',
    license='MIT',
    author='Qian Qin',
    author_email='qinqianhappy@gmail.com',
    description='DNase I seq QC pipeline',
    long_description = __doc__,
    scripts = ["gcap/GCAP", "gcap/pipeline-scripts/bed_duplicates.sh", "gcap/pipeline-scripts/script-tokenizer.py",
               "gcap/pipeline-scripts/runhotspot", "gcap/pipeline-scripts/macs2_spot.sh"],
    requires=["samflow", "jinja2", "argparse"],
    package_data = {"gcap" : ["static/*", "pipeline-scripts/*"]},
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'License :: OSI Approved :: Python Software Foundation License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules'
        ],
    )
