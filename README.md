# What is GCAP
> `G`lobal `C`hromatin `A`ccessibility `P`ipeline
> 
> X-nase(Dnase, Bnase, Cnase) Quality analysis pipeline
> 
---
## Prototype  Features

See the google docs [specification](https://docs.google.com/document/d/1dEcH3ezfODrL4ffMeqKU4YLBALBB18k61HAmCj9KozE/edit).[^1]

[^1]: ENCODE version

---
### Installation

Tested version of software bundled.

software    |version        | url |
:-----------| :-----------: | :-------|
python      | 3.3        | <https://www.python.org/ftp/python/3.3.3/Python-3.3.3.tgz> |
python      | 2.7        | <http://www.python.org/ftp/python/2.7.6/Python-2.7.6.tgz> | 
bwa         | 0.7.7        | <https://github.com/lh3/bwa> |
fastqStatsAndSubsample | 2 | 	<https://github.com/ENCODE-DCC/kentUtils> |
edwBamFilter | 2 | 	<https://github.com/ENCODE-DCC/kentUtils> |
samtools | 0.2.0 | https://github.com/samtools/samtools |
bedops | 2.4.2 | https://github.com/bedops/bedops/releases/download/v2.4.2/bedops_linux_x86_64-v2.4.2.tar.bz2 |
hotspot | 4 | https://github.com/rthurman/hotspot/archive/4.0.0.tar.gz |

bedtools |  2.17.0 | http://github.com/arq5x/bedtools/archive/v2.17.0.tar.gz |

Install `setuptools` before install GCAP. To avoid python version conflicts, users are recommended to use [virtualenv](https://raw.githubusercontent.com/pypa/virtualenv/1.9.X/virtualenv.py).

Install GCAP:

    git clone https://github.com/qinqian/GCAP
    python setup.py install

Outputs from GCAP
==========================
###QC metrics

###stored files

Reference
============
1. http://picard.sourceforge.net/explain-flags.html
2. [ENCODE 3 ChIP-seq pipeline](https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit)
3. https://github.com/ENCODE-DCC/uniformAnalysis
4. https://github.com/ENCODE-DCC/kentUtils

