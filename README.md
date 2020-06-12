# calpy

`calpy` is a CALPUFF model managing module in Python. 

[CALPUFF](http://www.src.com) is an advanced non-steady-state meteorological and air quality modeling system. Since its parameter file INP file (suffixed with ".inp") is old, ill-designed and hard to use, this python module provides functions to manage the INP file of CALPUFF model. 

If you want to manage other models or preprocessors INP files, you may need to implement them yourself. Since CALPUFF model's INP file is the most complicated one, other INP files management is easy.


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

The `data/` folder should contain your `calmet.dat` and options to run the CALPUFF model. Here I use a `data/winter/` folder to place all the data files and should set this prefix by InpGroup.set_data_dir("data/winter"). See the example for details.

Example:

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

calpy.run_calpuff(calgrp, disp=True)  # run CALPUFF

times, samples = calpy.load_samples(calgrp, timeZone="UTC+0800", filename="samples.csv")
timeseries = calpy.extract_receptors(calgrp)
simulates = timeseries.get_simulates(times)  # get simulated concentrations of each receptors

# do something
```

There are example csv files in this repository. 
