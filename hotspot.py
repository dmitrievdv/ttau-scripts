import numpy as np
import matplotlib.pyplot as plt
from fortmag import magnetosphere as mag
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.constants import year,pi,G
from scipy.optimize import fsolve

gridr = 4; gridtet = 5; levm = 5
Te=np.zeros((gridr,gridtet))
nh=np.zeros((gridr,gridtet)); ne=np.zeros((gridr,gridtet))
ni=np.zeros((gridr,gridtet,levm))
grid=np.zeros((gridr,gridtet,2))
T_hot=0.0

Rstar = 2.0
Mstar = 0.8
Tstar = 4000
rmi = 2.2
rmo = 3.0

grid = mag.init_grid(gridr, gridtet, rmi, rmo)

print(grid)

hotspot_dilution = mag.calc_hotspot_dilution(grid)

print(hotspot_dilution)
