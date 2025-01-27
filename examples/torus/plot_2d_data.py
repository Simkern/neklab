import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

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

def plot_2d_mesh(ax, pattern, only_edges=False, cwd='.'):

    x, y, vx, vy, vz, elmap, dt2d, metadata, nsteps = read_fields(pattern, only_mesh=True, cwd=cwd)

    nelf = metadata['nelf']
    # Loop to create the contours
    for i in range(nelf):
        xi = np.squeeze(x[:,:,i])
        yi = np.squeeze(y[:,:,i])
        plot_element(ax, xi, yi, only_edges)

def play_2d_data(pattern, cwd='.'):
    x, y, vx, vy, vz, elmap, dt2d, metadata, nsteps = read_fields(pattern, only_mesh=False, cwd=cwd)

    # Create the figure and axis
    fig, ax = plt.subplots()

    # Initialize the variables to store vmin and vmax for coloring
    vmin = vx[:,:,:,0].min()
    vmax = vx[:,:,:,0].max()
    dvmin = 1000
    dvmax = 0
    nelf = metadata['nelf']
    for i in range(nelf):
        dvmin = min(dvmin, (vx[:,:,i,:]).min())
        dvmax = max(dvmax, (vx[:,:,i,:]).max())

    # Define the update function for the animation
    for k in range(4):
        for j in range(0,nsteps,80):
            ax.clear()

            # Loop to create the contours
            for i in range(nelf):
                xi = np.squeeze(x[:,:,i])
                yi = np.squeeze(y[:,:,i])
                vxi = np.squeeze(vx[:,:,i,j])

                # Create the contour plot for each element
                contour = ax.contourf(xi, yi, vxi, vmin=dvmin, vmax=dvmax)

            # Add a color bar to the figure (only once)
            if j == 0 and k == 0:
                plt.colorbar(contour, ax= ax)
                ax.set_aspect('equal', 'box')

            ax.set_title(f'Plot {j + 1} of {nsteps}')  # Display current plot number (1-based index)

            plt.pause(0.001)

if __name__ == '__main__':
    play_2d_data('n2dtorus')
    