"""CALPUFF model managing module

Including:
- functions and class to read and write CALPUFF INP file
- functions to read CALPUFF LST file and draw plots
- functions for inversion

INP file is a parameters file for CALPUFF model. Only data within the
delimeters('!') are processed.

LST file is the CALPUFF output log file.

Author: Pavinberg
Email: pavin0702@gmail.com

"""

import sys
import re
import datetime
import pathlib
import subprocess
import itertools as it
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pygal


__all__ = [
    "InpGroups", "Series", "load_calpuff", "dump_calpuff", "run_calpuff",
    "get_coef", "extract_receptors", "extract_proportions",
    "load_samples", "draw_position", "draw_concentrations"
]


unitsVel = {
    "g/s": 1,
    "kg/h": 2,
    "lb/hr": 3,
    "tons/yr": 4,
    "m**3/s": 5,
    "m**3/min": 6,
    "mtons/yr": 7,
    "Bq/s": 8,
    "GBq/yr": 9
}

unitsConc = {
    "g/m**3": 1,
    "mg/m**3": 2,
    "ug/m**3": 3,
    "ng/m**3": 4,
    "odour": 5,
    "TBq/m**3": 6,
    "GBq/m**3": 7,
    "Bq/m**3": 8
}


class InpGroups:
    """A class to manage INP file information"""
    def __init__(self, idname="NULL", groups=None):
        self.id = idname  # module name of these groups
        self.groups = groups
        self._keys = {}  # key: variable - value: group number it belongs to
        # init keys from groups
        for i, group in enumerate(groups):
            for key in group.keys():
                if key in self._keys.keys():
                    raise ValueError(f"duplicate variable name {key} "
                                     f"in group{self._keys[key]} and {i}")
                self._keys[key] = i
        defaultParams = {
            # filenames -- Input Group 0
            "metFile": "calmet.dat",
            # time -- Input Group 1
            "startYear": 2019, "startMonth": 1, "startDay": 1, "startHour": 0,
            "endYear":   2019, "endMonth":   1, "endDay":   2, "endHour":   0,
            "timeZone": "UTC+0000",
            "timeStep": 3600,  # no larger than 3600
            # lambert projection -- Input Grroup 4
            "lambertlat1": "25N",
            "lambertlat2": "40N",
            "centerLat0": "31N",
            "centerLon0": "120E",
            # grids -- Input Grroup 4
            "nx": 60, "ny": 60, "nz": 11, "cellsize": 0.1,
            "zface": [0, 20, 40, 80, 160, 320, 640,
                      1000, 1500, 2000, 2500, 3000],  # length should be `nz`+1
            "xorigin": -3, "yorigin": -3
        }
        self.params = defaultParams
        self.dataDir = "data/"
        self.outDir = "output/"
        self.exeDir = "exe/"

    def __getitem__(self, varname):
        return self.groups[self._keys[varname]][varname]

    def __setitem__(self, varname, value):
        if isinstance(value, list):
            self.groups[self._keys[varname]][varname] = str(value)[1:-1]
        elif type(value) in [int, float, str]:
            self.groups[self._keys[varname]][varname] = value
        else:
            raise ValueError(
                f"\033[0;31;1mThere is an exception at line 114 in calpy.py"
                f"with type {type(value)}. Please contact the author"
                f"(pavin0702@gmail.com) with this information. "
                f"Thank You!\033[0m")

    def __len__(self):
        return len(self.groups)

    def set_data_dir(self, dirname: str):
        self.dataDir= dirname

    def set_output_dir(self, dirname: str):
        self.outDir = dirname

    def set_exe_dir(self, dirname: str):
        self.exeDir = dirname

    def set_params(self, params: dict = None):
        """Parse a dict to set some common parameters by alias names.
        Keys of params and their default values are:

        defaultParams = {
            # filenames -- Input Group 0
            "metFile": "calmet.dat",
            # time
            "startYear": 2019, "startMonth": 1, "startDay": 1, "startHour": 0,
            "endYear":   2019, "endMonth":   1, "endDay":   2, "endHour":   0,
            "timeZone": "UTC+0000",
            "timeStep": 3600,
            # lambert projection
            "lambertlat1": "25N",
            "lambertlat2": "40N",
            "centerLat0": "31N",
            "centerLon0": "120E",
            # grids
            "nx": 60, "ny": 60, "nz": 11, "cellsize": 0.1,
            "zface": [0, 20, 40, 80, 160, 320, 640,
                    1000, 1500, 2000, 2500, 3000],
            "xorigin": -3, "yorigin": -3
        }

        Args:

            params (dict/None): parameters to set.
                Default to None to set defaultParams.

        """
        if self.id != "CALPUFF":
            raise ValueError(f"this InpGroups id is {self.id}. You should "
                             f"load species only when the id is 'CALPUFF'.")
        if params:
            for key, value in params.items():
                if key in self.params.keys():
                    self.params[key] = value
                else:
                    raise KeyError(f"Unknown parameter \"{key}\"")
        realname = {
            # filenames -- Input Group 0
            "metFile": "METDAT",
            # time -- Input Group 1
            "startYear": "IBYR", "startMonth": "IBMO", "startDay": "IBDY", "startHour": "IBHR",
            "endYear":   "IEYR", "endMonth":   "IEMO", "endDay":   "IEDY", "endHour":   "IEHR",
            "timeZone": "ABTZ",
            "timeStep": "NSECDT",
            # lambert projection -- Input Group 4
            "lambertlat1": "XLAT1",
            "lambertlat2": "XLAT2",
            "centerLat0": "RLAT0",
            "centerLon0": "RLON0",
            # grids -- Input Grroup 4
            "nx": "NX", "ny": "NY", "nz": "NZ", "cellsize": "DGRIDKM",
            "zface": "ZFACE",
            "xorigin": "XORIGKM", "yorigin": "YORIGKM"
        }
        for key, value in self.params.items():
            real = realname[key]
            self.groups[self._keys[real]][real] = value
        self.groups[4]["IECOMP"] = self.params["nx"]
        self.groups[4]["IESAMP"] = self.params["nx"]
        self.groups[4]["JECOMP"] = self.params["ny"]
        self.groups[4]["JESAMP"] = self.params["ny"]
        self.groups[4]["ZFACE"] = str(self.params["zface"])[1:-1]
        self.groups[0]["METDAT"] = pathlib.Path(self.dataDir) \
            / self.params['metFile']
        outputFile = ["PUFLST", "CONDAT"]
        for param in outputFile:
            self.groups[0][param] = pathlib.Path(self.outDir) \
                / pathlib.Path(self.groups[0][param]).name

    def load_species_options(self, filename="species_options.csv",
                             unit="ug/m**3"):
        ''' load species basic information from a csv file.
        The header of csv file should be as follows:

        spec, modeled, emitted, deposited, droplets,
        conc_printed, conc_saved, dryflux_printed, dryflux_saved,
        wetflux_printed, wetflux_saved, massflux_saved

        '''
        if self.id != "CALPUFF":
            raise ValueError(f"this InpGroups id is {self.id}. You should "
                             f"load species only when the id is 'CALPUFF'.")

        filefullname = pathlib.Path(self.dataDir) / filename
        specOpt = pd.read_csv(filefullname, dtype=object)
        self.speciesList = specOpt[specOpt["%process"] == "1"]["spec"].to_numpy()
        self.speciesNum = len(self.speciesList)
        self.groups[1]["NSPEC"] = self.speciesNum
        self.groups[1]["NSE"] = self.speciesNum
        self.groups[5]["IPRTU"] = unitsConc[unit]

    def set_species(self, species):
        self.speciesList = np.array(species)
        self.speciesNum = len(self.speciesList)
        self.groups[1]["NSPEC"] = self.speciesNum
        self.groups[1]["NSE"] = self.speciesNum

    def load_species_decay(self, filename="species_decay.csv"):
        ''' load species decay half-life parameters.
        Use after load_species_options.
        '''
        if self.id != "CALPUFF":
            raise ValueError(f"this InpGroups id is {self.id}. You should "
                             f"load species only when the id is 'CALPUFF'.")
        if not hasattr(self, "speciesList"):
            raise ValueError(f"should call InpGroups.load_species_options before"
                             f" calling InpGroups.load_species_decay")
        colOrder = ['spec', 'half-life', 'mass-yield']
        filefullname = pathlib.Path(self.dataDir) / filename
        decay = pd.read_csv(filefullname, dtype=object)[colOrder]
        cond = [False] * decay.shape[0]
        speciesSet = set(self.speciesList)
        for i, row in decay.iterrows():
            if row['spec'] in speciesSet:
                cond[i] = True
        self.speciesDecay = decay[cond].to_numpy()
        self.groups[11]["NDECAY"] = self.speciesDecay.shape[0]
        self.groups[2]["MCHEM"] = 5

    def load_source_info(self, filename="source_info.csv"):
        if self.id != "CALPUFF":
            raise ValueError(f"this InpGroups id is {self.id}. You should "
                             f"load species only when the id is 'CALPUFF'.")
        colOrder = ['src', 'X', 'Y', 'stack-height', 'elevation',
                    'diameter', 'velocity', 'temperature', 'downwash']
        filefullname = pathlib.Path(self.dataDir) / filename
        srcInfo = pd.read_csv(filefullname, dtype=object)
        self.sourceInfo = srcInfo[srcInfo["%process"] == "1"][colOrder].to_numpy()
        self.sourceNum = self.sourceInfo.shape[0]
        self.groups[13]["NPT1"] = self.sourceNum

    def load_source_emission(self, filename="source_emission.csv", unit="g/s"):
        if self.id != "CALPUFF":
            raise ValueError(f"this InpGroups id is {self.id}. You should "
                             f"load species only when the id is 'CALPUFF'.")
        if not (hasattr(self, "speciesList") and hasattr(self, "sourceInfo")):
            raise ValueError(f"should call InpGroups.load_source_info and "
                             f"InpGroups.load_species_options before calling"
                             f"InpGroups.load_source_emission")
        filefullname = pathlib.Path(self.dataDir) / filename
        srcEmissAll = pd.read_csv(filefullname, dtype=object, index_col=0)
        try:
            srcEmiss = srcEmissAll.loc[self.sourceInfo[:, 0], self.speciesList]
        except KeyError:
            raise ValueError(
                "source emissions file does not contain all species that you specified."
            )
        self.sourceEmiss = srcEmiss.to_numpy().astype(float)
        self.groups[13]["IPTU"] = unitsVel[unit]

    def load_receptors_info(self, filename="receptors_info.csv"):
        """Read receptors information from a csv file

        Args:

            filename (str, optional): csv file name to read.
                Defaults to "receptors.csv".

        """
        colOrder = ["ID", "coorX", "coorY", "coorZ"]
        filefullname = pathlib.Path(self.dataDir) / filename
        self.receptors = pd.read_csv(filefullname, dtype=object)[colOrder].to_numpy()
        self.receptorsNum = self.receptors.shape[0]
        self.groups[21]["NREC"] = self.receptorsNum

    def load_area_source(self, filename="area_source.csv", unit="g/s"):
        ''' First line is basic information.
        Second line is polygon's X and third line is Y. 
        Forth line is emissions rate of the species
        '''
        if self.id != "CALPUFF":
            raise ValueError(f"this InpGroups id is {self.id}. You should "
                             f"load species only when the id is 'CALPUFF'.")
        filefullname = pathlib.Path(self.dataDir) / filename
        with open(filefullname, "r") as fp:
            self.areaSource = np.array(fp.readline().strip().split(','))
            self.areaSourceCoor = [None] * 2
            self.areaSourceCoor[0] = np.array(fp.readline().strip().split(','))
            self.areaSourceCoor[1] = np.array(fp.readline().strip().split(','))
            self.areaSourceEmiss = np.array(fp.readline().strip().split(','))
            if len(self.areaSourceEmiss) != self.speciesNum:
                raise ValueError(
                    f"area source file species number not corresponded to "
                    f"species specified "
                )
        self.groups[14]["NAR1"] = 1
        self.groups[14]["IARU"] = unitsVel[unit]

    def reset_source_emission(self, sourceEmiss):
        """Quick reset a new matrix of source emissions"""
        self.sourceEmiss = sourceEmiss

    def set_source_emission(self, sourceEmiss: np.ndarray, unit=None):
        if unit:
            self.groups[13]["IPTU"] = unitsVel[unit]
        if self.id != "CALPUFF":
            raise ValueError(f"this InpGroups id is {self.id}. You should "
                             f"load species only when the id is 'CALPUFF'.")
        self.sourceEmiss = sourceEmiss

    def set_area_source_emission(self, areaEmiss):
        self.areaSourceEmiss = areaEmiss

    def dump_source_emission(self, filename):
        with open(filename, "w") as fp:
            fp.write(f"src,{','.join(self.speciesList)}\n")
            srcEmiss = []
            for line in self.sourceEmiss:
                srcEmiss.append(list(map(lambda x: f"{x:.4f}", line)))
            mat = np.hstack((self.sourceInfo[:, 0:1], srcEmiss))
            for line in mat:
                fp.write(",".join(line) + '\n')

    def printGroup(self, groupnum: int):
        if groupnum > len(self.groups):
            raise ValueError(
                "group number out of range. There are {len(self.groups)} groups."
            )
        for key, value in self.groups[groupnum].items():
            if key[:3] != "END":
                print(f"{key:<8} = {value}")

    def dump(self, fp, groupnum: int = -1):
        ''' output parameters of group: `groupname` to file pointer fp.
        If groupnum == -1: output all groups.
        '''
        if groupnum > len(self.groups):
            raise ValueError(
                "group number out of range. There are {len(self.groups)} groups."
            )
        if groupnum == -1:
            # dump all groups
            for group in self.groups:
                for key, value in group.items():
                    fp.write(f"! {key} = {value} !\n")
                fp.write("!END!\n")
        else:
            for key, value in self.groups[groupnum].items():
                fp.write(f"! {key} = {value} !\n")
            fp.write("!END!\n")

    def raw_set(self, **kwargs):
        ''' set parameters in groups with real name '''
        for key, value in kwargs.items():
            if key in self._keys.keys():
                self.groups[self._keys[key]][key] = value
            else:
                raise ValueError(f"variable name {key} not found in InpGroups")

    def get_time_info(self):
        """get start date, end date and time step

        Returns:

            datetime.datetime: start date
            datetime.datetime: end date
            datetime.timedelta: time step
        """
        startDate = datetime.datetime(
            self.params["startYear"], self.params["startMonth"],
            self.params["startDay"], self.params["startHour"]
        )
        endDate = datetime.datetime(
            self.params["endYear"], self.params["endMonth"],
            self.params["endDay"], self.params["endHour"]
        )
        timeStep = datetime.timedelta(seconds=self.params["timeStep"])
        return startDate, endDate, timeStep


def basic_read(content: str) -> OrderedDict:
    ''' read parameters from a template string '''
    items = re.findall("!.*?!", content)
    params = OrderedDict()
    for item in items:
        pair = re.split("!|=", item[1:-1])
        pair = list(map(str.strip, pair))
        if pair[0] == "END":
            ...
            # global ENDCount
            # params[f"END{ENDCount}"] = ENDCount
            # ENDCount += 1
        else:
            params[pair[0]] = pair[1]
    return params


def load_calpuff(filename="utils/calpuff_template.inp") -> InpGroups:
    """load calpuff_template.inp as template

    Args:
        filename (str, optional): calpuff template INP file.

    Returns:
        calpy.InpGroups: a class to manage INP file parameters.
    """
    SPECIAL_GROUPS = [3, 7, 8]
    with open(filename, "r") as fp:
        contents = fp.read().split("INPUT GROUP")[1:]
        if len(contents) != 22:
            print("Should have 22 groups. Check the file input.")
            sys.exit()
        grps = [None] * 22
        for i, content in enumerate(contents):
            if i in SPECIAL_GROUPS:
                grps[i] = OrderedDict()
            else:
                grps[i] = basic_read(content)
        return InpGroups(idname="CALPUFF", groups=grps)


def cut_line(line, maxlength=80):
    pos = [i for i, ch in enumerate(line) if ch == ","]
    nx = []
    old = 0
    for idx, nxt in zip(pos, np.roll(pos, shift=-1)):
        if nxt - old >= maxlength:
            nx.append(line[old:idx+1])
            old = idx + 1
    nx.append(line[old:])
    return "\n".join(nx)


def dump_calpuff(groups: InpGroups, filename="calpuff.inp"):
    ''' output a calpuff.inp file'''
    if groups.id != "CALPUFF":
        raise ValueError("groups' id should be 'CALPUFF'")

    def group3dump(fp, i):
        if hasattr(groups, "speciesList"):
            for spec in groups.speciesList:
                fp.write(f"! CSPEC = {spec} !  !END!\n")
            for spec in groups.speciesList:
                fp.write(f"! {spec} = 1, 1, 0, 0, 0 !\n")
            fp.write("!END!\n")

    def group5dump(fp, i):
        if hasattr(groups, "speciesList"):
            for spec in groups.speciesList:
                fp.write(f"! {spec} = 1, 1, 0, 0, 0, 0, 0 !\n")
        groups.dump(fp, 5)

    def group7dump(fp, i):
        fp.write("!END!\n")

    def group8dump(fp, i):
        fp.write("!END!\n")

    def group10dump(fp, i):
        fp.write("!END!\n")
        groups.dump(fp, 10)

    def group11dump(fp, i):
        groups.dump(fp, 11)
        if hasattr(groups, "speciesDecay"):
            for line in groups.speciesDecay:
                fp.write(f"! {line[0]} = {line[1]}, {line[2]} ! !END!\n")
                #fp.write("!END!\n")
        fp.write("!END!\n")

    def group13dump(fp, i):
        groups.dump(fp, 13)
        if hasattr(groups, "sourceInfo"):
            srcEmiss = groups.sourceEmiss
            for line, emiss in zip(groups.sourceInfo, srcEmiss):
                emissStr = list(map(lambda x: f"{x:.4f}", emiss))
                xline = ','.join(np.concatenate([line[1:], emissStr]))
                xline = f"! X = {xline} !\n"
                if len(xline) > 80:
                    xline = cut_line(xline)

                fp.write(
                    f"! SRCNAM = {line[0]} !\n"
                    f"{xline}"
                    f"! ZPLTFM = 0 !\n"
                    f"! FMFAC = 1.0 ! "
                    f"!END!\n")

    def group14dump(fp, i):
        groups.dump(fp, 14)
        if hasattr(groups, "areaSource"):
            fp.write(f"! SRCNAM = {groups.areaSource[0]} !\n")
            xline = (f"! X = {','.join(groups.areaSource[1:])},"
                     f" {','.join(groups.areaSourceEmiss)} ! !END!\n")
            # if len(xline) > 80:
            #         xline = cut_line(xline)
            fp.write(xline)
            fp.write(f"! SRCNAM = {groups.areaSource[0]} !\n")
            fp.write(f"! XVERT = {','.join(groups.areaSourceCoor[0])} !\n")
            fp.write(f"! YVERT = {','.join(groups.areaSourceCoor[1])} ! !END!\n")

    def group21dump(fp, i):
        groups.dump(fp, 21)
        if hasattr(groups, "receptors"):
            for line in groups.receptors:
                fp.write(f"! X = {line[1]}, {line[2]}, {line[3]} ! !END!\n")

    def default(fp, i):
        groups.dump(fp, i)

    groups_func = {
        3: group3dump,
        5: group5dump,
        7: group7dump,
        8: group8dump,
        10: group10dump,
        11: group11dump,
        13: group13dump,
        14: group14dump,
        21: group21dump
    }

    # !!CAREFUL!! the first line of title is very important. 
    # It must contain 3 parts with lengths of 16+16+(48) repectively.
    # Use spaces for pacing. Make sure the 32 < length <= 80.
    title = "CALPUFF.INP     7.01             Created by inp-module\n"
    #        |---------------|---------------|----------------------------------------------|
    #        0              16              32                                             80
    
    # filefullname = pathlib.Path(groups.exeDir) / filename
    with open(filename, "w") as fp:
        fp.write(title+"\n\n\n")
        for i in range(len(groups)):
            groups_func.get(i, default)(fp, i)


def run_calpuff(groups: InpGroups, exe="calpuff.exe", disp=True):
    """Run CALPUFF model with the groups

    Args:

        groups (InpGroups): InpGroups with parameters to run the model

        exe (str, optional): calpuff executable file name. 
            Defaults to "exe/calpuff.exe".

    Returns:

        bool: True if run successfully.
    """
    dump_calpuff(groups)
    exefilename = pathlib.Path(groups.exeDir) / exe
    ret = subprocess.run(str(exefilename), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if disp:
        print(ret.stdout.decode("utf8"))
        print(ret.stderr.decode("utf8"))
    if b"HALTED" in ret.stderr or b"ERROR" in ret.stderr or b"severe" in ret.stderr:
        print("Error when running calpuff")
        sys.exit()
    return True


def get_coef(groups: InpGroups, simulates: np.ndarray):
    """get coefficient which is (C/Q)_sim

    Args:

        groups (calpy.InpGroups): paramenters

        simulates (numpy.ndarray): simulated result created by the model.
            dimension: [species x receptors x (source+1)]

    Returns:

        numpy.ndarray: [species x receptors x source]

    """
    simu = simulates[:, :, :-1]  # last vectors are total, should be excluded
    src = groups.sourceEmiss.transpose()
    # for every species
    mp = map(lambda sim, emiss: np.nan_to_num(sim / emiss), simu, src)
    return np.array(list(mp))


class Series:
    """storage time series as a 4-dimensions numpy array.

    period_number x source_number x receptors_number x species_number
    """
    def __init__(self, timeSeries: np.ndarray, startDate: datetime.datetime,
                 timeStep: datetime.timedelta):
        self.timeSeries = timeSeries
        self.startDate = startDate
        self.timeStep = timeStep
        self.shape = timeSeries.shape
        self.sourceNum = self.shape[1] - 1
        self.speciesNum = self.shape[3]
        self.iteri = 0

    def __iter__(self):
        self.iteri = 0
        return self

    def __next__(self):
        if self.iteri < self.speciesNum:
            frame = self.timeSeries[self.iteri]
            self.iteri += 1
            return frame
        raise StopIteration

    def get_simulates(self, times: list):
        # species x receptors x sources
        simulates = np.empty((self.speciesNum, len(times), self.sourceNum + 1))
        for recep, time in enumerate(times):
            simulates[:, recep, :] = self.timeSeries[time, :, recep, :].transpose()
        return simulates


def extract_receptors(groups: InpGroups, filename=None):
    """read calpuff.lst file and get concentrations 
    information for each receptors

    """
    if groups.id != "CALPUFF":
        raise ValueError("groups' id should be 'CALPUFF'")
    if not filename:
        filename = groups["PUFLST"]
    startDate, endDate, timeStep = groups.get_time_info()
    timeFramesNum = (endDate - startDate) // timeStep  # should be divided exactly
    sourceNum = groups.sourceNum
    if hasattr(groups, "areaSource"):
        sourceNum += 1
    recepNum = groups.receptorsNum
    speciesNum = groups.speciesNum
    seriesLen = timeFramesNum * (sourceNum + 1)

    series = np.empty((seriesLen, recepNum, speciesNum))
    with open(filename, "r") as fp:
        blocks = it.filterfalse(
            lambda x: x[0],
            it.groupby(fp, lambda x: x.startswith(" Receptor No."))
        )
        blocks = it.islice(blocks, 1, None)  # drop the first one
        for iframe, block in zip(range(seriesLen), blocks):
            for irecp, line in zip(range(recepNum), block[1]):
                record = line.split()[1:]
                series[iframe, irecp] = np.array(list(map(float, record)))
    series = series.reshape(timeFramesNum, sourceNum + 1, recepNum, speciesNum)
    # series = series.transpose(3, 0, 1, 2)  # species x time x source x receptors
    return Series(series, startDate, timeStep)


def extract_proportions(groups: InpGroups):
    """ Get species proportions of each source 
    and the total source emission rates

    Returns:

    numpy.ndarray: proportions, [sources x species]
    numpy.ndarray: total emission of every source [sources x 1]
    """
    srcEmiss = groups.sourceEmiss
    total = np.array([np.sum(src) for src in srcEmiss]).reshape(-1, 1)
    return total, np.nan_to_num(srcEmiss / total)


def load_samples(groups: InpGroups, filename="samples.csv", timeZone="UTC+0800"):
    if not re.match(r"UTC[+|-](0\d|1[012])", timeZone):
        raise ValueError(f"unknown format of UTC time zone '{timeZone}''. "
                         f"Should be like 'UTC+0800' ")
    filefullname = pathlib.Path(groups.dataDir) / filename
    df = pd.read_csv(filefullname, index_col=0)
    try:
        concentrations = df.loc[groups.receptors[:, 0], groups.speciesList].to_numpy()
    except KeyError:
        raise ValueError(
            "samples file does not contain all species that you specified."
        )

    def getDatetime(dateStr: str):
        dateList = np.array(dateStr.split("-")).astype(int)
        return datetime.datetime(*dateList)

    def getIndex(date: datetime.datetime):
        startDate, _, timeStep = groups.get_time_info()
        return (date - startDate) // timeStep

    timeDiff = int(timeZone[3:6]) - int(groups.params['timeZone'][3:6])
    dates = df['time'].apply(getDatetime)
    times = np.array(list(map(getIndex, dates)))
    times = times - timeDiff  # index of time corresponded to receptors
    return times, concentrations


def draw_position(groups, filename="position.svg", dir="images",
                  drawSource=True, drawReceptors=True):
    ''' plot sources' and receptors' position '''
    from pygal import Config
    config = Config()
    config.title = "Position"
    config.show_lengend = False
    config.height = 500
    config.stroke = False
    config.dots_size = 2.5
    plot = pygal.XY(
        config, x_title="X/m", y_title="Y/m", show_x_guides=True,
        print_labels=True, style=pygal.style.styles["default"]
        (value_font_size=2)
    )
    if drawSource:
        src = groups.sourceInfo
        meta = []
        for line in src:
            name = line[0]
            coor = line[1:3].astype(float) * 1000
            meta.append({'value': tuple(coor), 'label': name})
        plot.add("sources", meta)
    if drawReceptors:
        recep = groups.receptors
        meta = []
        for line in recep:
            name = line[0]
            coor = line[1:3].astype(float) * 1000
            meta.append({'value': tuple(coor), 'label': name, 'dots_size': 2})
        plot.add("receptors", meta)
    plot.render_to_file(pathlib.Path(dir) / filename)


def draw_concentrations(groups: InpGroups, species=None,
                        filenmae="images/concentration.png",
                        inputfile=None, drawSource=True, drawReceptors=True):
    ''' read the calpost output about concentrations and plot a heatmap'''
    xorigin, yorigin = groups.params['xorigin'], groups.params['yorigin']
    nx, ny = groups.params['nx'], groups.params['ny']
    cellsize = groups.params['cellsize']
    minx, miny = 0, 0
    maxx, maxy = nx, ny

    if drawSource:
        coors = groups.sourceInfo[:, 1:3].astype(float).T
        srcx = (coors[0] - xorigin) / cellsize
        srcy = (coors[1] - yorigin) / cellsize
        minx, miny = np.min(srcx), np.min(srcy)
        maxx, maxy = np.max(srcx), np.max(srcy)
    if drawReceptors:
        coors = groups.receptors[:, 1:3].astype(float).T
        recpx = (coors[0] - xorigin) / cellsize
        recpy = (coors[1] - yorigin) / cellsize
        minx, miny = min(minx, np.min(recpx)), min(miny, np.min(recpy))
        maxx, maxy = max(maxx, np.max(recpx)), max(maxy, np.max(recpy))

    margin = int(max(nx, ny) * 0.05)
    minx, miny = max(int(minx) - margin, 0), max(int(miny) - margin, 0)
    maxx, maxy = min(int(maxx) + margin, nx), min(int(maxy) + margin, ny)

    if inputfile is None:
        if species is None:
            raise ValueError("parse species name as argument")
        species = species.lower()
        concfile = f"output/calpost/rank(0)_{species}_12hr_conc.csv"
    else:
        concfile = inputfile
    with open(concfile) as fp:
        for i in range(6):
            fp.readline()
        df = pd.read_csv(fp, header=None)
        conc = df.iloc[:, 2].to_numpy()
    concgrid = conc.reshape(nx, ny)
    plt.pcolor(concgrid, cmap=plt.cm.Blues)
    plt.colorbar()
    if drawSource:
        plt.scatter(srcx, srcy, label="sources", c="red", alpha=0.5, s=10)
        plt.legend()
    if drawReceptors:
        plt.scatter(recpx, recpy, label="receptors", c="orange", alpha=0.5, s=10)
        plt.legend()
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)
    plt.title(species)
    plt.show()
