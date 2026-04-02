import numpy as np
import physical_geodesy as pg

Nmax = 1000
lat = 10.0

P = pg.physical_geodesy.fcm_afl_fuku_py(Nmax, lat)

for m in range(Nmax + 1):
    for n in range(m, Nmax + 1):
        print(f"P({n:3d},{m:3d}) = {P[n,m]:15.5E}")
