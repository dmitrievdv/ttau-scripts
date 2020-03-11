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

model = 'hart94_90_70_22-30'  # -- Name of the model
suffix = '_hym'

have_observ = True # -- If True will get observation data from file ./observation.dat
					# -- This file should include two columns of data -- first for velocites
					# -- and second for normalized intensity
try: # -- Checking if file ./observation.dat exist by trying to open it
	observ = open('observation.dat', 'r')
	observ_data = np.loadtxt(observ)
	observ_vels = observ_data[:,0]
	observ_prof = observ_data[:,1]
	print(observ_data.shape)
	observ.close()
except:
	have_observ = False # -- If it doesn't exist set have_observ to False


try: 
    population_parameters = rm.read_parameters(model) # -- Trying to read model parameters
except rm.NoSuchModelError as err: 
    print("Can't find data for model "+err.model_name)
    population_parameters = None # -- If failed set parameters to none
        
print(population_parameters) # -- Printing parameters for checking
field_type = population_parameters['field_type'] # -- Setting field type

u,l = 3,2 # -- Setting upper and lower level of the line

popul_model = population_parameters['populations'] + suffix # -- Setting model name for populations
try: # -- Trying to read the populations
    (interp_grid, Te_atgrid, nh_atgrid, 
     ne_atgrid, n_u_atgrid, n_l_atgrid) = rm.read_populations_file(popul_model, u, l, field_type)
except rm.NoSuchModelError as err: # -- Self-explanatory
    print("Can't find populations for model "+err.model_name)
except rm.NoDataForLevel as err: # -- Self-explanatory
    print("Level "+str(err.level)+" to big for model "+err.model_name)
except rm.UpLevelLowerThenLowLevel: # -- Self-explanatory
    print('Up level is lower than low level')

n,m,mz = 100,100,100 # -- Setting accuracy:
					 # -- 1. n  - number of frequencies (velocities) on which the profile 
					 # --         will be calculated
					 # -- 2. m  - Number of points across one dimension of a 2d grid on the picture plane.
					 # --         Total number of grid points on the picture plane: m^2
					 # -- 3. mz - Number of points on z-axis (line of sight integration)

i,psi,alpha = 15,0,0 # -- Orientation parameters:
					 # -- 1. i     - angle between star rotational axis and line of sight
					 # -- 2. alpha - angle between star rotational axis and the field axis
					 # -- 3. psi   - phase angle of rotation of the field axis around rotational axis

Rstar = population_parameters['Rstar'] # -- Setting radius of the star
Mstar = population_parameters['Mstar'] # -- Setting mass of the star
Tstar = population_parameters['Tstar'] # -- Setting temperature of the star

first_border = population_parameters['first_border'] # -- Setting first border of the field
second_border = population_parameters['second_border'] # -- Setting second border of the field
in_cut = population_parameters['in_cut'] # -- Setting inner field-cutting radius
out_cut = population_parameters['out_cut'] # -- Setting outer field-cutting radius

Mdot = population_parameters['Mdot'] # -- Setting acretion rate

vrot = 15 # -- Setting rotational velocity of the star on the equator

# a = np.array([1,1]) # -- What?! I don't remember why is it here and is it even needed...

prof.init_star(Rstar, Mstar, Tstar, vrot) # -- Star initialization
prof.init_field(field_type, Mdot, first_border, second_border, in_cut, out_cut) # -- Field initialization
# prof.init_custom_line(10830e-8, 5, 3, 4, 7.28591e-13, 1e4) # -- Custom line initialization
														   # -- arguments: (wave length, upper stat. weight, lower stat. weight
														   # --    atom mass, absorbtion coeficient, temperature for the absorbtion coef.)
prof.init_hydrogen_line(u, l) # -- Hydrogen line initialization
frequencies, nu_0 = prof.init_frequencies(n, 1) # -- Frequency grid initialization. Second argument controls borders.
												# -- Examples: 1 - [-vesc:vesc], 0.5 - [-0.5*vesc:0.5*vesc]

prof.init_orientation(i, psi, alpha) # -- Init field orientation
field_borders, grid, dS = prof.init_field_borders(m, grid_type = 'polar') # -- Calculation of the line of sight field borders 
                                                                          # -- for each picture plain grid point

# -- Initialization of populations interpolation: start
init_state = {'dipole' : prof.init_dipole_state, 'cone' : prof.init_cone_state} 
init_state[field_type](n_u_atgrid, n_l_atgrid, Te_atgrid, ne_atgrid,
     nh_atgrid, interp_grid)
# -- Initialization of populations interpolation: end

# -- Interface creation: start
fig, axes = plt.subplots(1,2) 
plt.subplots_adjust(right = 0.89, bottom=0.25) 
fig.suptitle('i = ' + str(i) + ', m = ' +  str(m)+', model: '+model)
	

axes[1].axis('equal')
axes[1].set_ylim(-out_cut, out_cut)
axes[1].set_xlim(-out_cut*1.33333, out_cut*1.333333)
axes[1].set_adjustable("box")

axcolor = plt.axes([0.9, 0.3, 0.02, 0.5])

# axS = plt.axes([0.2, 0.17, 0.7, 0.03], facecolor='gray')
axBS = plt.axes([0.2, 0.1, 0.15, 0.07], facecolor='gray')
axS =  plt.axes([0.35, 0.1, 0.2, 0.07], facecolor='gray')
axB = plt.axes([0.55, 0.1, 0.35, 0.07], facecolor='gray')
axI = plt.axes([0.6, 0.03, 0.1, 0.07], facecolor='gray')
axBI = plt.axes([0.7, 0.03, 0.2, 0.07], facecolor='gray')
axM = plt.axes([0.2, 0.03, 0.1, 0.07], facecolor='gray')
axBM = plt.axes([0.3, 0.03, 0.2, 0.07], facecolor='gray')
# slS =  Slider(axS, r'$\nu, km/s$', -500, 500, valinit=0)
But = Button(axB, 'Plot map!!!', color = '0.85', hovercolor='1')
Buts = Button(axBS, 'SaveProf', color = '0.85', hovercolor='1')
Txts = TextBox(axS, 'file')
Buti = Button(axBI, 'Change i', color ='0.85', hovercolor='1')
Txti = TextBox(axI, 'i = ')
Butm = Button(axBM, 'Change m', color ='0.85', hovercolor='1')
Txtm = TextBox(axM, 'm = ')

# -- Interface creation: start


# -- Unneccesary solving of radiation transfer equation for one grid point: start 
dz, z_coord, r, n_u, n_l, ut, stark, v_z, dv_d, alphaline, k_lu, tau, S_ul = prof.check_point_freq(0, 0, nu_0, mz)

output = open('z.dat', 'w')
	
for ii in range(100):
	output.write(str(z_coord[i])+' '+str(dz[i])+ ' '+ str(r[i])+ ' '+ str(n_u[i])+ ' '+ str(n_l[i])+
	        ' '+ str(ut[i])+ ' '+ str(stark[i])+ ' '+ str(v_z[i])+ ' '+ str(dv_d[i])+ ' '+ str(alphaline[i])+
	        ' '+ str(k_lu[i])+ ' '+ str(tau[i])+ ' '+ str(S_ul[i])+ '\n')

output.close()
# -- Unneccesary solving of radiation transfer equation for one grid point: end

stark = True
profile, emission_map, emission = prof.calc_profile(frequencies, field_borders, grid, dS, mz, no_stark = not stark) # -- Profile calculation



freq_vel = 0
fake_mouse_event = MouseEvent('fake', fig.canvas, 0,0) # -- Faking click on central frequency



def recalculate(): # -- Recalculate the profile 
	global field_borders, grid, dS, profile, emission_map, emission
	# prof.init_orientation(i, psi, alpha)
	field_borders, grid, dS = prof.init_field_borders(m, grid_type = 'polar')
	profile, emission_map, emission = prof.calc_profile(frequencies, field_borders, grid, dS, mz, no_stark = not stark)
	fig.suptitle('i = ' + str(i) + ', m = ' +  str(m)+', model: '+model)

def update_m(val): # -- Update m (picture plane grid accuracy)
	global m, profile, grid, dS, field_borders, emission_map
	try:
		m = int(Txtm.text)
	except:
		msg.showerror(title="Input error", message="Wrong m input")
		return
	# prof.init_orientation(i, psi, alpha)
	# prof.init_field(field_type, Mdot, first_border, second_border, in_cut, out_cut)
	recalculate()
	# print(profile)
	redraw(1)
	redraw_prof(fake_mouse_event)

def update_i(val): # -- Update declination angle
	global i, field_borders, grid, dS, profile, emission_map
	try:
		i = int(Txti.text)
	except:
		msg.showerror(title="Input error", message="Wrong i input")
		return
	# fig.suptitle('Wait....')
	# fig.canvas.draw_idle()
	prof.init_orientation(i, psi, alpha)
	# prof.init_field(field_type, Mdot, first_border, second_border, in_cut, out_cut)
	recalculate()
	# print(profile)
	redraw(1)
	redraw_prof(fake_mouse_event)

def redraw_prof(event): # -- Redrawing the profile
	global freq_vel
	try:
		if(event.inaxes in [axes[0]]):
			freq_vel = event.xdata
	except:
		pass
		return
	axes[0].clear()
		
	# fig.canvas.draw_idle()
	freq = nu_0 - freq_vel/3e5*nu_0
	closest_ind = int(np.abs(frequencies - freq).argmin())
	closest_freq = frequencies[closest_ind]
	if(freq >= closest_freq):
		freq_interval = (frequencies[closest_ind], 
							frequencies[closest_ind+1])
		prof_interval = (profile[closest_ind], 
							profile[closest_ind+1])
	else:
		freq_interval = (frequencies[closest_ind-1], 
								frequencies[closest_ind])
		prof_interval = (profile[closest_ind-1], 
							profile[closest_ind])

	delta_freq = freq_interval[1] - freq_interval[0]
	delta_prof = prof_interval[1] - prof_interval[0]
	prof_val = (freq - freq_interval[0])/delta_freq*delta_prof + prof_interval[0]

	if have_observ:
		axes[0].plot(observ_vels, observ_prof, 'k--', linewidth = 0.5)
	axes[0].plot(-(frequencies-nu_0)/nu_0*3e5, profile, 'b-')
	axes[0].plot([-(freq-nu_0)/nu_0*3e5, -(freq-nu_0)/nu_0*3e5],
						 [1, prof_val], 'r-')
	fig.canvas.draw_idle()
	fig.canvas.set_window_title(str(u)+'->'+str(l))

def redraw(val): #  -- Redrawing the emission map
	global emission_map
	axes[1].clear()
	axcolor.clear()
	axes[1].axis('equal')
	axes[1].set_ylim(-out_cut, out_cut)
	axes[1].set_xlim(-out_cut*1.33333, out_cut*1.333333)
	axes[1].set_adjustable("box")

	freq = nu_0 - freq_vel/3e5*nu_0
	dummy, emission_map, emission = prof.calc_profile([freq], field_borders, grid, dS, mz, no_stark = not stark)
	axes[1].plot(np.cos(np.linspace(0, 2*np.pi, 100)), np.sin(np.linspace(0, 2*np.pi, 100)), 'r--')#, zorder = 1)
	cmap = plt.cm.Greys
	em_map_plot = axes[1].pcolormesh(grid[1,:,:], grid[0,:,:], np.ma.log10(emission_map[:-1,:-1]), cmap = cmap)#, zorder = 2)
	cbar = fig.colorbar(em_map_plot, cax = axcolor, ax = axes[1], format='%.1f')
	fig.canvas.draw_idle()


def save(val): # Saving the profile in text file
	filename = Txts.text
	filename = filename.replace(' ', '')
	if filename == '':
		filename = model + suffix +'_i'+str(i)+'_m'+str(m)+'_vrot'+str(int(vrot))
	output = open('profiles/'+filename+'.dat', 'w')
	for freq, prof, emm in zip(frequencies, profile, emission):
		output.write(str(-(freq-nu_0)/nu_0*3e5)+' '+str(prof)+' '+str(emm)+' '+str(freq)+'\n')

	output.close()

	map_output = open('maps/'+filename+'{:+d}.dat'.format(int(freq_vel)), 'w')
	# freq = nu_0 - freq_vel/3e5*nu_0
	map_output.write('# '+str(freq_vel) + '\n')
	for k in range(m):
		for j in range(m):
			# print(str(grid[1,j,i]), str(grid[0,j,i]), str(emission_map[j,i]))
			map_output.write('{:>7.3f} {:>7.3f} {:>11.4e}\n'.format(grid[1,j,k], grid[0,j,k], emission_map[j,k]))
			# map_output.write(str(grid[1,j,i])+' '+str(grid[0,j,i])+' '+str(emission_map[j,i])+'\n')
		map_output.write('\n')

	map_output.close()

redraw_prof(fake_mouse_event)
redraw(1)
fig.canvas.mpl_connect("button_press_event", redraw_prof)
# slS.on_changed(redraw_prof)
But.on_clicked(redraw)
Buti.on_clicked(update_i)
Butm.on_clicked(update_m)
Buts.on_clicked(save)
plt.show()

