import numpy as np
import readmodel as rm
import matplotlib
import tkinter.messagebox as msg
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt 
from matplotlib.widgets import Slider, TextBox, Button
from matplotlib.backend_bases import MouseEvent
from fortprof import profile as prof
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.constants import year,pi,G
from scipy.optimize import fsolve
from matplotlib.colors import ListedColormap
import os

models = ['UXOR_70_78_12-17', 'UXOR_70_78_15-25', 'UXOR_70_93_12-17', 'UXOR_70_93_15-25', 
          'UXOR_75_78_12-17', 'UXOR_75_78_15-25', 'UXOR_75_93_12-17', 'UXOR_75_93_15-25', 
          'UXOR_80_78_12-17', 'UXOR_80_78_15-25', 'UXOR_80_93_12-17', 'UXOR_80_93_15-25']
profdir = 'UXOR'
suffixes = ['']
angles = [70]
lines = [(3,2), (4,2), (5,2)]

try:
    os.mkdir(profdir)
except: 
    pass

for line in lines:
    u,l = line
    for mainmodel in models:
        for suffix in suffixes:
            for i in angles:
                model = mainmodel + suffix
                print(model)
                population_parameters = rm.read_parameters(model)
                field_type = population_parameters['field_type']
                popul_model = population_parameters['populations']
                (interp_grid, Te_atgrid, nh_atgrid,
                    ne_atgrid, n_u_atgrid, n_l_atgrid) = rm.read_populations_file(popul_model, u, l, field_type)
                n,m,mz = 100,100,100
                psi,alpha = 0,0
                vrot = 150

                Rstar = population_parameters['Rstar']
                Mstar = population_parameters['Mstar']
                Tstar = population_parameters['Tstar']

                Mdot = population_parameters['Mdot']
                # if(abs(np.log10(Mdotmod/Mdot)) > 1e-2): 
                #     print('Mdot error')
                #     quit() 

                first_border = population_parameters['first_border']
                second_border = population_parameters['second_border']
                in_cut = population_parameters['in_cut']
                out_cut = population_parameters['out_cut']

                prof.init_star(Rstar, Mstar, Tstar, vrot)
                prof.init_field(field_type, Mdot, first_border, second_border, in_cut, out_cut)
                # prof.init_custom_line(10830e-8, 5, 3, 4, 7.28591e-13, 1e4)
                prof.init_hydrogen_line(u, l)
                frequencies, nu_0 = prof.init_frequencies(n, 1)

                prof.init_orientation(i, psi, alpha)
                field_borders, grid, dS = prof.init_field_borders(m, grid_type = 'polar')
                # plt.pcolormesh(grid[1,:,:], grid[0,:,:], field_borders[:-1,:-1]), cmap = 'Greys')

                # prof.init_orientation(i,psi,alpha)
                # field_borders, grid, dS = prof.init_field_borders(m)

                init_state = {'dipole' : prof.init_dipole_state, 'cone' : prof.init_cone_state}
                init_state[field_type](n_u_atgrid, n_l_atgrid, Te_atgrid, ne_atgrid,
                                            nh_atgrid, interp_grid)

                profile, emission_map, emission = prof.calc_profile(frequencies, field_borders, grid, dS, mz, no_stark = False)

                # plt.plot(-(frequencies-nu_0)/nu_0*3e5, profile, 'b-')
                # plt.show()

                output = open('{:s}/{:s}_l{:d}{:d}_i{:d}_prof.dat'.format(profdir, model, u, l, i), 'w')
                print(profile[0], frequencies[0])
                for vel, profil in zip(frequencies, profile):
                    vel = -(vel-nu_0)/nu_0*3e5
                    output.write('{:15f} {:15f}\n'.format(vel, profil))
                output.close()
                # quit()