import numpy as np
import traceback
import matplotlib.pyplot as plt
import readmodel as readm
from fortmag import magnetosphere as mag
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.constants import year,pi,G
from scipy.optimize import fsolve
import os


model_data_directory = 'models/data/'
model_popul_directory = 'models/popul/'

os.makedirs(model_data_directory, exist_ok = True)
os.makedirs(model_popul_directory, exist_ok = True)

gridr = 3; gridtet = 20; levm = 15



Te=np.zeros((gridr,gridtet))
nh=np.zeros((gridr,gridtet)); ne=np.zeros((gridr,gridtet))
ni=np.zeros((gridr,gridtet,levm))
ni_loc = ni
grid=np.zeros((gridr,gridtet,2))
T_hot=0.0



def init_d_grid(rmi, rmo):
    fig, axes = plt.subplots(2,1)
    plt.subplots_adjust(right = 0.89, bottom=0.25) 
    x = [0,0.2,1]; y = [0.25,0.5,0.25]
    step = np.array(list(zip(x,y)))
    def draw_grid(step):
        axes[0].clear()
        axes[1].clear()
        axes[1].axis('equal')
        axes[1].set_xlim(0, rmo*1.1)
        axes[0].set_ylim(0.025,1)
        axes[0].set_xlim(-0.1,1.1)
        axes[1].set_adjustable("box")
        d_grid = mag.calc_d_grid(gridtet, step.transpose())
        print(d_grid)
        grid = mag.init_grid(gridr, rmi, rmo, mtet=gridtet, grid_d = d_grid)
        print(grid[:,0,1])
        cart_x = grid[:,:,0]*np.sin(grid[:,:,1])**3
        cart_y = grid[:,:,0]*np.sin(grid[:,:,1])**2*np.cos(grid[:,:,1])
        axes[0].plot(step.transpose()[0,:], step.transpose()[1,:], 'b-')
        for plot_x, plot_y in zip(cart_x, cart_y):
            axes[1].plot(plot_x,plot_y,'k--',lw=0.5)
        axes[1].plot(cart_x, cart_y, 'ko', ms = 2)
        fig.canvas.draw_idle()
        return d_grid
    def redraw_grid(event):
        try:
            if(event.inaxes in [axes[0]]):
                if(event.button == 1):
                    if(event.xdata < 0):
                        y[0] = event.ydata
                    elif(event.xdata > 1):
                        y[1] = event.ydata
                    else:
                        x.append(event.xdata)
                        y.append(event.ydata)
                if(event.button == 3):
                    min_ind = np.argmin(abs(x-event.xdata)[2:])
                    x.pop(min_ind+2)
                    y.pop(min_ind+2)
                step = np.array(list(zip(x,y)))
                step = step[step[:,0].argsort()]
                d_grid = draw_grid(step)
        except:
            traceback.print_exc()
            return
    print('kek')
    d_grid = draw_grid(step)

    fig.canvas.mpl_connect("button_press_event", redraw_grid)
    # plt.show()
    step = np.array(list(zip(x,y)))
    step = step[step[:,0].argsort()]
    d_grid = draw_grid(step)
    # print(d_grid)
    return d_grid

# model_names = [#'RZPsc_9_85_6-10', 'RZPsc_10_90_8-12', 'RZPsc_9_85_8-10', 'RZPsc_9_85_10-12',
#                'RZPsc_9_90_8-12', 'RZPsc_9_90_8-10', 'RZPsc_9_90_10-12',
#                'RZPsc_9_90_6-10', 'RZPsc_9_85_8-12', 'RZPsc_10_85_8-10',
#                'RZPsc_10_85_10-12', 'RZPsc_10_85_6-10', 'RZPsc_10_85_8-12',
#                'RZPsc_10_90_8-10', 'RZPsc_10_90_10-12', 'RZPsc_10_90_6-10']

# model_names = ['UXOR_80_80_12-17', 'UXOR_75_80_15-25', 'UXOR_70_80_15-25',
#                'UXOR_80_80_12-17', 'UXOR_75_80_12-17', 'UXOR_70_80_12-17']#, 'hart94_80_70_22-30', 'hart94_90_70_22-30']

model_names = ['hart94_90_75_22-30', 'hart94_85_75_22-30', 'hart94_80_75_22-30', 
               'hart94_90_75_42-50', 'hart94_85_75_42-50', 'hart94_80_75_42-50']
               # 'hart94_80_75_22-40', 'hart94_85_75_22-40', 'hart94_90_75_22-40',
               # 'hart94_80_75_42-50', 'hart94_85_75_42-50', 'hart94_90_75_42-50']
# model_names = ['UXOR_70_78_12-17', 'UXOR_70_78_15-25', 'UXOR_70_93_12-17', 'UXOR_70_93_15-25', 
# 			   'UXOR_75_78_12-17', 'UXOR_75_78_15-25', 'UXOR_75_93_12-17', 'UXOR_75_93_15-25', 
# 			   'UXOR_80_78_12-17', 'UXOR_80_78_15-25', 'UXOR_80_93_12-17', 'UXOR_80_93_15-25']
# 			   #'hart94_80_70_22-30', 'hart94_82_70_22-30', 'hart94_85_70_22-30', 'hart94_87_70_22-30', 'hart94_90_70_22-30', 
               


# ['PhiLeo_10_75_20-25', 'PhiLeo_10_80_20-25', 'PhiLeo_95_75_20-25',
#                'PhiLeo_95_80_20-25', 'PhiLeo_9_75_20-25', 'PhiLeo_9_80_20-25',
#                'PhiLeo_10_75_15-20', 'PhiLeo_10_80_15-20', 'PhiLeo_95_75_15-20',
#                'PhiLeo_95_80_15-20', 'PhiLeo_9_75_15-20', 'PhiLeo_9_80_15-20']


#hart94_85_8K
# model_name = input('Model name: ')

for model_name in model_names:
    try:
        model_parameters = readm.read_parameters(model_name)
    except readm.NoSuchModelError as err:
        print("Can't find model "+err.model_name)
        continue

    Te=np.zeros((gridr,gridtet))
    nh=np.zeros((gridr,gridtet)); ne=np.zeros((gridr,gridtet))
    ni=np.zeros((gridr,gridtet,levm))
    ni_loc = ni
    ni_rev = ni
    # grid=np.zeros((gridr,gridtet,2))

    Rstar, Mstar, Tstar = (model_parameters[key] for key in ['Rstar', 'Mstar', 'Tstar'])
    T_norm, Mdot = (model_parameters[key] for key in ['Tmax', 'Mdot'])
    T_norm = 6000
    rmi, rmo = (model_parameters[key] for key in ['first_border', 'second_border'])
    # Mdot = 1e-7
    # T_hot = (Mdot*2e30/3e7*6.67e-11*Mstar*2e30/Rstar/7e8*(1-2/(rmi+rmo))/4/pi/(Rstar*7e8)**2/0.1/5.67e-8)**0.25
    # print('T_hot:', T_hot)
    
    # print(Rstar, Mstar, Tstar)
    # print(T_norm, Mdot)
    # print(rmi, rmo)

    print("Model: ", model_name, " T_max: ", T_norm, " lg(Mdot): ", np.log10(Mdot), 'Mag. size: '+str(rmi)+'-'+str(rmo))

    # continue
    # Mdot = 10**(Mdot)

    d_grid = init_d_grid(rmi,rmo)
    grid = mag.init_grid(gridr, rmi, rmo, mtet=gridtet, grid_d = d_grid)
    print(d_grid)
    print("Grid created")
    print(grid.shape)
    print(grid[0,:,1])

    cart_x = grid[:,:,0]*np.sin(grid[:,:,1])**3
    cart_y = grid[:,:,0]*np.sin(grid[:,:,1])**2*np.cos(grid[:,:,1])

    # for plot_x, plot_y in zip(cart_x, cart_y):
    #     plt.plot(plot_x,plot_y,'k--')

    # for plot_x, plot_y in zip(cart_x.T, cart_y.T):
    #     plt.plot(plot_x,plot_y,'k--')

    # plt.plot(cart_x, cart_y, 'ko', ms = 5)
    # plt.axis('scaled')
    # plt.show()

    # continue
    T = np.full((gridr,gridtet), T_norm)

    # def iterate(T, f, slow = 0):
       
    #     nh, ne, ni = mag.calc_populations(grid, 0.1, Mdot, Mstar, Rstar, Tstar, T_hot, 7, loc=True)
    #     f_new = ne/nh
    #     # print(np.amax(np.abs(1-f_new/f)))
    #     if(slow != 0):
    #         ne = (ne - f*nh)/slow + f*nh
    #     T = get_T(grid, nh, ne, T_norm, Rstar)
    #     return T, ne/nh

    # ionization = np.full((gridr, gridtet), 1)


    # T, ionization = iterate(T, ionization)
    
    r = grid[:,:,0]*np.sin(grid[:,:,1])**2
    hot_spot, hot_T, hot_d = (model_parameters[key] for key in ['hot_spot', 'Thot', 'Dhot'])
    # hot_T = 3000
    # hot_d = 0.1
    print('Hot spot: ', hot_spot)
    if(hot_spot):
        hot = hot_T*np.exp(-(r-1)/hot_d)
        nh, ne, ni = mag.calc_populations(grid, 0.2, Mdot, Mstar, Rstar, Tstar, T_hot, levm, model=model_name)#, loc=False)
    else:
        
        nh, ne, ni = mag.calc_populations(grid, T_norm, Mdot, Mstar, Rstar, Tstar, T_hot, levm, model=model_name,
                                             loc = True, use_hymera = False, nh_cooling = True, model_dir=model_popul_directory,
                                             lc_ionization = True)
        os.rename(model_popul_directory+'/'+model_name+'_popul.dat', 
                  model_popul_directory+'/'+model_name+'_nhcool-stat_loc_popul.dat')

        nh, ne, ni = mag.calc_populations(grid, T_norm, Mdot, Mstar, Rstar, Tstar, T_hot, levm, model=model_name,
                                             loc = True, use_hymera = True, nh_cooling = False, model_dir=model_popul_directory,
                                             lc_ionization = True)
        os.rename(model_popul_directory+'/'+model_name+'_popul.dat', 
                  model_popul_directory+'/'+model_name+'_nonstat_loc_popul.dat')
        os.rename(model_popul_directory+'/'+model_name+'_loc_popul.dat', 
                  model_popul_directory+'/'+model_name+'_stat_loc_popul.dat')

        

    
    # for l in range(7):
    #     plt.title(str(l+1)+' level')    
    #     for plot_x, n_l in zip(cart_x, ni[:,:,l]/ni_loc[:,:,l]):
    #         plt.plot(plot_x, n_l, 'k-')
    #     plt.show()

    mag.clean_memory()

    # quit()



