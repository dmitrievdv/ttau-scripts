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

u,l = 3,2

models = ['hart94_75_70_22-30', 'hart94_77_70_22-30', 'hart94_80_70_22-30','hart94_82_70_22-30','hart94_85_70_22-30','hart94_90_70_22-30']
Mdots = [10**-7.5, 10**-7.7, 10**-8.0, 10**-8.2, 10**-8.5, 10**-9]
suffixes = ['', '_simpcool', '_simpcool_stat']
angles = [15, 45, 75]
eqwidth = np.zeros((len(suffixes), len(angles), len(Mdots)))

eqsuf = 0
for suffix in suffixes:
    eqa = 0
    filename = 'eqwidth'+suffix
    output = open('profiles/'+filename+'.dat', 'w')
    for i in angles:
        eqmd = 0
        for mainmodel,Mdot in zip(models, Mdots):
            model = mainmodel + suffix
            print(model)
            population_parameters = rm.read_parameters(model)
            field_type = population_parameters['field_type']
            popul_model = population_parameters['populations']
            (interp_grid, Te_atgrid, nh_atgrid,
                ne_atgrid, n_u_atgrid, n_l_atgrid) = rm.read_populations_file(popul_model, u, l, field_type)
            n,m,mz = 100,100,100
            psi,alpha = 0,0
            vrot = 15

            Rstar = population_parameters['Rstar']
            Mstar = population_parameters['Mstar']
            Tstar = population_parameters['Tstar']

            Mdotmod = population_parameters['Mdot']
            if(abs(np.log10(Mdotmod/Mdot)) > 1e-2): 
                print('Mdot error')
                quit() 

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
            # plt.plot(-(frequencies-nu_0)/nu_0*3e5, np.ma.masked_where(profile > 1, profile), 'b-')
            # plt.show()

            deltav = abs(frequencies[1]-frequencies[0])/nu_0*3e5
            cureqwidth = np.ma.sum(np.ma.masked_where(profile < 1, profile) - 1)*deltav
            eqwidth[eqsuf, eqa, eqmd] = cureqwidth
            print('eqwidth: '+str(cureqwidth))
            output.write(str(Mdot)+' '+str(cureqwidth)+'\n')
            eqmd += 1 
        output.write('\n') 
        # plt.plot(np.log10(Mdots), np.log10(eqwidth[eqsuf, eqa, :]), 'bo')
        # plt.show()
        eqa += 1
    output.close()
    eqsuf += 1

