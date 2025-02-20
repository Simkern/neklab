import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from read_2d_data import read_fields

def plot_element(ax, x, y, only_edges=False, color='k'):
    l = x.shape[0]
    lm = x.shape[0] - 1
    if not only_edges:
        for i in range(lm):
            ax.plot(x[i,:l], y[i,:l], c=color)
            ax.plot(x[:l,i], y[:l,i], c=color)
        color = 'r'
    else:
        color = color
    for i in [0, lm]:
        ax.plot(x[i,:l], y[i,:l], c=color)
        ax.plot(x[:l,i], y[:l,i], c=color)
    ax.set_aspect('equal', 'box')

def plot_2d_mesh(ax, x, y, only_edges=False, color='k'):

    nelf = x.shape[-1]
    # Loop to create the contours
    for i in range(nelf):
        xi = np.squeeze(x[:,:,i])
        yi = np.squeeze(y[:,:,i])
        plot_element(ax, xi, yi, only_edges, color=color)
    ax.set_aspect('equal', 'box')

def plot_2d_orientation(ax, x, y, offset=2, draw_elnum=True):

    nelf = x.shape[-1]
    xa = np.sum(x, axis=(0,1))/np.prod(x.shape[:2])
    ya = np.sum(y, axis=(0,1))/np.prod(x.shape[:2])
    # Loop to create the contours
    o = offset
    for i in range(nelf):
        xi = np.squeeze(x[:,:,i])
        yi = np.squeeze(y[:,:,i])
        plot_element(ax, xi, yi, only_edges=True)
        ax.plot(xi[o:-o,o], yi[o:-o,o],c='red',linewidth=2)
        if draw_elnum:
            ax.text(xa[i], ya[i], f'{i+1}', color='black', fontsize=12)
    ax.scatter(x[o,o,:],y[o,o,:],s=30,c='red',marker='o')
    ax.set_aspect('equal', 'box')

def plot_2d_fld(ax, x, y, fld, istep=0, draw_elements=False, draw_mesh=False):

    nelf = x.shape[-1]
    if not draw_mesh:
        only_edges=True

    vmin = fld[:,:,:,istep].min()
    vmax = fld[:,:,:,istep].max()

    for i in range(nelf):
        xi = np.squeeze(x[:,:,i])
        yi = np.squeeze(y[:,:,i])
        vxi = np.squeeze(fld[:,:,i,istep])
        c = ax.contourf(xi, yi, vxi, vmin=vmin, vmax=vmax)
        if draw_elements:
            plot_element(ax, xi, yi, only_edges)
    ax.set_aspect('equal', 'box')
    fig = plt.gcf()
    c.set_clim(vmin, vmax)
    cbar = fig.colorbar(c, ax=ax)

def animate_2d_fld(ax, data, step=80, if_half=False):

    # Initialize the variables to store vmin and vmax for coloring
    nelf = data.x.shape[-1]
    nsteps = data.vx.shape[-1]
    dvmin = data.vx.min()
    dvmax = data.vx.max()
    T = sum(data.dt)

    # Define the update function for the animation
    for k in range(4):
        for j in range(0,nsteps,step):
            ax.clear()

            # Loop to create the contours
            for i in range(nelf):
                xi = np.squeeze(data.x[:,:,i])
                yi = np.squeeze(data.y[:,:,i])
                vxi = np.squeeze(data.vx[:,:,i,j])

                # Create the contour plot for each element
                contour = ax.contourf(xi, yi, vxi, vmin=dvmin, vmax=dvmax)
                if if_half:
                    contour = ax.contourf(-xi, yi, vxi, vmin=dvmin, vmax=dvmax)

            # Add a color bar to the figure (only once)
            if j == 0 and k == 0:
                ax.set_aspect('equal', 'box')

            t = sum(data.dt[:j])
            tp = t/T*100
            ax.set_title(f't = {t:5.3f}     ({tp:3.0f} % T)\n', fontsize=20)  # Display current plot number (1-based index)
            plt.axis('off')
            plt.pause(0.001)

def play_file(pattern, cwd='.'):
    data, meta, nsteps = read_fields(pattern, only_mesh=False, cwd=cwd)

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(15, 15))

    animate_2d_fld(ax, data, if_half = meta.if_half)

if __name__ == '__main__':
    #play_file('n2dtorus')
    #play_file('prod/run/Wo_040.0/Q_0.180/n2dtorus')
    #play_file('prod/run/Wo_025.0/Q_0.350/n2dtorus')
    play_file('prod/run/Wo_035.0/Q_0.280/n2dtorus')