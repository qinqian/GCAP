from configparser import ConfigParser, NoSectionError, NoOptionError
import os
import time


class NoTreatmentData(Exception):
    pass

class Conf(object):
    def __init__(self, conf, args):
        self._verbose_level = 1
        self._conf = ConfigParser()

        if not os.path.exists(conf):
            raise IOError("No such config file: %s" % repr(conf))

        self._conf.read(conf)


        self.threads = args.threads
        self.pe = args.pe

        self.input = args.input
        self.target_dir = args.out
        self.species = args.species
        
        self.id = args.name
        
        self._conf.set("basics", "time", time.strftime("%Y-%m-%d"))
        self._conf.set("basics", "species", self.species)
        self._conf.set("basics", "id", self.id)
        self._conf.set("basics", "input", self.input)
        self._conf.set("basics", "output", os.path.abspath(self.target_dir))

        if args.species in ["hg19", "hg38"]:
            self._conf.set("macs2", "species", "hs")
        if args.species in ["mm9", "mm10"]:
            self._conf.set("macs2", "species", "mm")

        f = open('%s.conf'%(self.id), 'w')
        self._conf.write(f)
        self.root_dir = os.path.dirname('%s.conf'%(self.id))
        f.close()


    def set_option(self, verbose_level=1):
        """
        :type verbose_level: int
        verbose_level:  1 - only show fatal errors          (quiet mode)
                        2 - show fatal errors and warnings  (normal mode)
                        3 - show workflow details           (verbose mode)
                        4 - debug mode                      (debug mode)
        """
        self._verbose_level = verbose_level

    def get(self, section, option, default = None):
        try:
            return self._conf.get(section, option)
        except NoOptionError:
            if default:
                return default
            else:
                raise

    def get_path(self, section, option):
        return self.to_abs_path(self.get(section, option))

    def items(self, section):
        try:
            return self._conf.items(section)
        except NoSectionError:
            if self._verbose_level >= 2:
                print("Warning: No such section: ", section)
                print("This will return a empty dict")
            return {}

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self,value):
        self._id = value

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, value):
        self._threads = value

    @property
    def pe(self):
        return self._pe

    @pe.setter
    def pe(self, value):
        '''setting Pair End state, True for PE, False for SE
        '''
        self._pe = value

    @property
    def target_dir(self):
        return self._target_dir

    @target_dir.setter
    def target_dir(self, value):
        self._target_dir = value

    @property
    def prefix(self):
        return os.path.join(self.target_dir, self.id)

    @property
    def json_prefix(self):
        return os.path.join(self.category("json"), self.id)

    @property
    def latex_prefix(self):
        return os.path.join(self.category("latex"), self.id)

    @property
    def hotspot_merge_final_prefix(self):
        return os.path.join(self.prefix + "_merge_all-final", self.id)

    @property
    def hotspot_reps_final_prefix(self):
        return [ os.path.join(f + "-final", self.id + "_treat_rep" + str(i+1)) for i, f in enumerate(self.treatment_targets) ]

    @property
    def hotspot_reps_final_5M_prefix(self):
        return [ os.path.join(f + "_5M_sort-final", self.id + "_treat_rep" + str(i+1) + "_5M_sort") for i, f in enumerate(self.treatment_targets) ]

    @property
    def hotspot_starch_input(self):
        """ bed.starch should have a different directory from output """
        return [ os.path.join(self.category("hotspot_starch_input"), self.id + "_treat_rep" + str(i+1)) for i, f in enumerate(self.treatment_targets) ]

    @property
    def hotspot_merge_starch(self):
        """ bed.starch should have a different directory from output """
        return os.path.join(self.category("hotspot_starch_input"), self.id)

    @property
    def treatment_pairs(self):
        """
        one to one in single end mode,  original, target
        two to one in pair end mode [original_pair1, original_pair2], target
        """
        return list(zip(self.treatment_raws, self.treatment_targets))

    @property
    def treatment_pairs_pe(self):
        return list(zip(self.treatment_raws, self.treatment_pair_targets["pairs"]))

    def to_abs_path(self, path):
        abs_path = path
        if not os.path.isabs(path):
            abs_path = os.path.join(self.root_dir, abs_path)
        return abs_path

    @property
    def treatment_bases(self):
        return [os.path.basename(i) for i in self.treatment_targets]

    @property
    def treatment_raws(self):
        """
        single end data separate by , for replicates
        pair end data separate by ; for replicates , for pairs
        """
        if self.get("basics", "input"):
            if not self.pe:
                return [self.to_abs_path(i.strip()) for i in self.get("basics", "input").split(",")]
            else:
                data_list = []
                for i in self.get("basics", "input").split(";"):
                    data_list.append([ self.to_abs_path(j.strip()) for j in i.split(",") ])
                return data_list
        else:
            raise NoTreatmentData

    @property
    def treatment_targets(self):
        if not self.pe:
            return self.treatment_single_targets
        else:
            return self.treatment_pair_targets["reps"]

    @property
    def treatment_pair_data(self):
        return self.treatment_pair_targets["pairs"]

    @property
    def treatment_single_targets(self):
        return [os.path.join(self.target_dir,
            self.id + "_treat_rep" + str(num+1)
        ) for num in range(len(self.treatment_raws))]

    @property
    def treatment_pair_targets(self):
        '''pairs: for [[rep1_pair1, rep1_pair2]], 
        usually for evaluating read quality
        reps: for [rep1, rep2],
        usually for mapping pair end data
        '''
        return {"pairs": [ [os.path.join(self.target_dir, self.id + "_treat_rep" + str(num+1)) + "pair1",
                  os.path.join(self.target_dir, self.id + "_treat_rep" + str(num+1)) + "pair2"]
                 for num in range(len(self.treatment_raws)) ],
                "reps": [os.path.join(self.target_dir,
                    self.id + "_treat_rep" + str(num+1)
                ) for num in range(len(self.treatment_raws))]}

    def category(self, category_name):
        target_path = os.path.join(self.target_dir, category_name)
        if not os.path.exists(target_path):
            os.makedirs(target_path)
        return target_path


    @property
    def peakcalltool(self):
        return self.get("tool", "peak_calling").strip().lower()
