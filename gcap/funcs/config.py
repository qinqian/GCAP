from configparser import ConfigParser, NoSectionError, NoOptionError
import os


class NoTreatmentData(Exception):
    pass

class Conf(object):
    def __init__(self, conf):
        self._verbose_level = 1
        self._conf = ConfigParser()
        if not os.path.exists(conf):
            raise IOError("No such config file: %s" % repr(conf))
        self._conf.read(conf)
        self.root_dir = os.path.dirname(conf)

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
        return self.get("Basis", "id")

    @property
    def target_dir(self):
        return self.get("Basis", "output")

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
    def treatment_bam(self):
        """
        inital input BAMs or SAMs, BEDs, separated by comma, no matter PE or SE
        """
        return [self.to_abs_path(i.strip()) for i in self.get("Basis", "treat").split(",")]

    @property
    def treatment_bed(self):
        """
        initial input reads BAMs, SAMs, BEDs, separated by comma, no matter PE or SE
        """
        return [self.to_abs_path(i.strip()) for i in self.get("Basis", "treat").split(",")]

    @property
    def treatment_raws(self):
        """
        single end data separate by , for replicates
        pair end data separate by ; for replicates , for pairs
        """
        if self.get("Basis", "treat").strip():
            if self.seq_type == "se":
                return [self.to_abs_path(i.strip()) for i in self.get("Basis", "treat").split(",")]
            elif self.seq_type == "pe":
                data_list = []
                for i in self.get("Basis", "treat").split(";"):
                    data_list.append([ self.to_abs_path(j.strip()) for j in i.split(",") ])
                return data_list
            elif self.seq_type.startswith("bam") or self.seq_type.startswith("sam") or self.seq_type.startswith("bed"):
                return self.treatment_bam
        else:
            raise NoTreatmentData

    @property
    def treatment_targets(self):
        if self.seq_type == "se":
            return self.treatment_single_targets
        elif self.seq_type == "pe":
            return self.treatment_pair_targets["reps"]
        elif self.seq_type.startswith("bam") or self.seq_type.startswith("sam") or self.seq_type.startswith("bed"):
            return self.treatment_single_targets

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
    def maptool(self):
        return self.get("tool", "mapping").strip().lower()

    @property
    def peakcalltool(self):
        return self.get("tool", "peak_calling").strip().lower()

    @property
    def seq_type(self):
        return self.get("Basis", "sequence_type").strip().lower()