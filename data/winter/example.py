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
