"""
GCAP
-----

GCAP is a python module for wrapping different tools to quality control
DNase-seq and other open chromatin sequencing data.

GCAP is easy to start
`````````````````````

.. code:: python

   GCAP -i "test_pair2.fastq,test_pair1.fastq;test_pair2.fastq,test_pair1.fastq" -n test -o test -s hg19 --pe --threads 18 


We suggest using virtualenv to avoid python version library conflicts.

Links
`````

* `website <https://github.com/qinqian/GCAP>`_

"""

from setuptools import setup
from configparser import SafeConfigParser
import os

version = str(1.1)

def generate_defaults():
    """Tries to generate the default run configuration files from the
    paths specified in chilin.conf.  The idea is that an admin, who is
    installing ChiLin system-wide could define the defaults which will allow
    users to avoid generating/editing conf files!!
    """

    cf = SafeConfigParser()
    cf.add_section('basics')
    ls = ['id', 'time', 'species', 'input', 'output']
    for fld in ls:
        cf.set('basics', fld, '$'+fld.upper())
    #SET the chilin version number
    cf.set('basics', "version", version)

    #read in the chilin.conf file--look first for developer defaults
    if os.path.exists('gcap.conf.filled'):
        cf.read('gcap.conf.filled')
    else:
        cf.read('gcap.conf')

    #write the template file!
    f = open(os.path.join('gcap','static','gcap.conf.filled'),'w')
    cf.write(f)
    f.close()


def main():
    if os.path.exists("gcap.conf"):
        generate_defaults()

    setup(
        name='GCAP',
        version=version,
        packages=['gcap', "gcap.funcs", "samflow"],
        url='https://github.com/qinqian/GCAP',
        license='GNU',
        author='Qian Qin',
        author_email='qianqind@gmail.com',
        description='DNase I seq QC pipeline',
        long_description = __doc__,
        scripts = ["gcap/GCAP",
                   "gcap/glue/dac_se_read_quality",
                   "gcap/glue/dac_pe_read_quality",
                   "gcap/glue/eap_run_bwa_se",
                   "gcap/glue/eap_run_bwa_pe",
                   "gcap/glue/dac_bam_se_post_filter",
                   "gcap/glue/dac_bam_pe_post_filter",
                   "gcap/glue/eap_dnase_stats",
                   "gcap/glue/eap_run_hotspot",
                   "gcap/glue/hotspot.py",
                   "gcap/glue/dac_pbc",
                   "gcap/glue/bdg2bw",
                   "gcap/glue/dac_spot",
                   "gcap/glue/dac_macs2_spot.sh",
                   "gcap/glue/eap_narrowPeak_to_bigBed",
                   "gcap/glue/eap_broadPeak_to_bigBed",
                  ],
        package_data = {"gcap" : ["static/*", "glue/*as"]},
        install_requires=['jinja2','argparse'],
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

if __name__ == "__main__":
    main()
