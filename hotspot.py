import numpy as np
import matplotlib.pyplot as plt
from fortmag import magnetosphere as mag
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.constants import year,pi,G
from scipy.optimize import fsolve

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


gridr = 6; gridtet = 40; levm = 5
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

d_grid = init_d_grid(rmi,rmo)
grid = mag.init_grid(gridr, rmi, rmo, mtet=gridtet, grid_d = d_grid)

print(grid)

mean_temp = mag.calc_mean_surf_temp(4000, 10000, grid)

cart_x = grid[:,:,0]*np.sin(grid[:,:,1])**3
cart_y = grid[:,:,0]*np.sin(grid[:,:,1])**2*np.cos(grid[:,:,1])
cart_r = (cart_x**2 + cart_y**2)**0.5

fig, axes = plt.subplots(1,1)

for plot_r, T in zip(cart_r, mean_temp):
    # plt.plot(plot_r, star_dilut, 'k--')
    plt.plot(plot_r, T, 'k-')
    # plt.axhline(4000, color = 'k', linestyle = '--')
    # plt.axhline(10000, color = 'k', linestyle = '--')

plt.show()
