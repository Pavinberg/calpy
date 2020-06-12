# calpy -- CALPUFF model managing module in Python

[CALPUFF](http://www.src.com) is an advanced non-steady-state meteorological and air quality modeling system. Since its parameter file INP file (suffixed with ".inp") is old, ill-designed and hard to use, this python module provides functions to manage the INP file of CALPUFF model. This module needs python >= 3.6.

The main idea is to use a class called `InpGroup` to operate the parameters by its methods. In the example below you will see operations like setting directory path, setting species information, setting source information and so on. In fact these operations just save the information into the data structure contained in the InpGroup class. When you call `calpy.run_calpuff` function, it will output all the parameters to a file called `calpuff.inp` and run the `calpuff.exe` in the `exe/` directory.

The information is set in the csv files in `data/` directory. As long as you set them right, `calpy` will do the things left, e.g. setting number of species and sources. If you need to do multiple tasks, you can create a sub-directory in `data/` just like the example below, where I use a `data/winter/` directory to place all the data files and set this prefix by InpGroup.set_data_dir("data/winter"). See code below for details.

Before using the module, you need to place your `calpuff.exe` into the `exe/` directory and `calmet.dat` into the `data/winter/` directory.

The module is too simple for practical purposes. So you may need to understand this module before using it. It only support Lambert projection and use a base `calpuff_template.inp` file in the `utils/` directory to set default parameters value. So you may need to alter to your own `calpuff_template.inp` file.

If you want to manage other preprocessors INP files, you may need to implement them yourself. Since CALPUFF model's INP file is the most complicated one, other INP files management is easy.

The project directory should be like:

```
.
├── calpy.py
├── clean.py
├── data
│   └── winter
│       ├── area_source.csv
│       ├── calmet.dat
│       ├── receptors_info.csv
│       ├── samples.csv
│       ├── source_emission.csv
│       ├── source_info.csv
│       ├── species_decay.csv
│       └── species_options.csv
├── exe
│   └── calpuff.exe
├── images
├── output
│   └── winter
├── utils
│   └── calpuff_template.inp
└── yourCode.py
```

Example codes:

```python
import calpy

params = {
    # CALMET file
    "metFile": "calmet.dat",
    # time -- Input Group 1
    "startYear": 2019, "startMonth":  1, "startDay":  23, "startHour":  0,
    "endYear":   2019, "endMonth":    1, "endDay":    23, "endHour":   12,
    "timeZone": "UTC+0000",
    "timeStep": 3600,  # no larger than 3600
    # lambert projection -- Input Grroup 4
    "lambertlat1": "25N",
    "lambertlat2": "40N",
    "centerLat0": "31N",
    "centerLon0": "120E",
    # grids -- Input Grroup 4
    "nx": 60, "ny": 60, "nz": 11, "cellsize": 1,
    "zface": [0, 20, 40, 80, 160, 320, 640,
		      1000, 1500, 2000, 2500, 3000],  # length should be `nz`+1
    "xorigin": -30, "yorigin": -30
}
group = calpy.load_calpuff()
group.set_data_dir("data/winter")
group.set_output_dir("output/winter")
group.set_exe_dir("exe/")
group.set_params(params)
group.load_species_options(unit="ug/m**3")
print("species:\n", group.speciesList)
group.load_species_decay()
# print("decay:\n", group.speciesDecay)
group.load_source_info()
# print("source:\n", group.sourceInfo)
group.load_source_emission(unit="g/s")
# print("emissn:\n", group.sourceEmiss)
group.load_receptors_info()
# print("receptors:\n", group.receptors)
#group.load_area_source()

calpy.run_calpuff(group, disp=True)  # run CALPUFF

times, samples = calpy.load_samples(group, timeZone="UTC+0800", filename="samples.csv")
timeseries = calpy.extract_receptors(group)
simulates = timeseries.get_simulates(times)  # get simulated concentrations of each receptors

# do something
```

This will load all the information set in `data/winter/` and run the model. All the output of CALPUFF model will be put into `output/winter/` directory.

There are example csv files in this repository. The first columns of `source_info.csv` and `species_options.csv` are "%process" which indicate whether to add the source or species into the INP file, 1 for yes and 0 for no.

Area source's csv file is too simple and crude for I didn't have time to re-design it. You may want to rewrite this part on your own.

`clean.py` is used to clean the annoying files created by CALPUFF.
