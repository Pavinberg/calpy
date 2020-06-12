import os

removeList = ["luse.clr", "pftrak.dat", "qagrid.bna", "qaluse.grd",
              "qapnts.dat", "qarecd.dat", "qarecg.dat", "qaterr.grd"]

for file in removeList:
    os.remove(file)

print("Cleaned")