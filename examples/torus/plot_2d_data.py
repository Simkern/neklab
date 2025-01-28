import numpy as np
import matplotlib.pyplot as plt

from read_2d_data import read_fields

def plot_element(ax, x, y, only_edges=False):
    l = x.shape[0]
    lm = x.shape[0] - 1
    if not only_edges:
        for i in range(lm):
            ax.plot(x[i,:l], y[i,:l], c='k')
            ax.plot(x[:l,i], y[:l,i], c='k')
        c = 'r'
    else:
        c = 'k'
    for i in [0, lm]:
        ax.plot(x[i,:l], y[i,:l], c=c)
        ax.plot(x[:l,i], y[:l,i], c=c)
    ax.set_aspect('equal', 'box')

def plot_2d_mesh(ax, x, y, only_edges=False):

    nelf = x.shape[-1]
    # Loop to create the contours
    for i in range(nelf):
        xi = np.squeeze(x[:,:,i])
        yi = np.squeeze(y[:,:,i])
        plot_element(ax, xi, yi, only_edges)
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

def animate_2d_fld(ax, x, y, fld, step=80):

    # Initialize the variables to store vmin and vmax for coloring
    nelf = x.shape[-1]
    nsteps = fld.shape[-1]
    dvmin = fld.min()
    dvmax = fld.max()

    # Define the update function for the animation
    for k in range(4):
        for j in range(0,nsteps,step):
            ax.clear()

            # Loop to create the contours
            for i in range(nelf):
                xi = np.squeeze(x[:,:,i])
                yi = np.squeeze(y[:,:,i])
                vxi = np.squeeze(fld[:,:,i,j])

                # Create the contour plot for each element
                contour = ax.contourf(xi, yi, vxi, vmin=dvmin, vmax=dvmax)

            # Add a color bar to the figure (only once)
            if j == 0 and k == 0:
                plt.colorbar(contour, ax=ax)
                ax.set_aspect('equal', 'box')

            ax.set_title(f'Plot {j + 1} of {nsteps}')  # Display current plot number (1-based index)

            plt.pause(0.001)

def play_file(pattern, cwd='.'):
    x, y, vx, vy, vz, elmap, dt2d, metadata, nsteps = read_fields(pattern, only_mesh=False, cwd=cwd)

    # Create the figure and axis
    fig, ax = plt.subplots()

    animate_2d_fld(ax, x, y, vx)

if __name__ == '__main__':
    play_file('n2dtorus')