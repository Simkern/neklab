import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

from read_2d_data import read_binary_file

vx_list = []
vy_list = []
vz_list = []

file_paths = [ 'n2dtorus001.fld', 'n2dtorus002.fld', 'n2dtorus003.fld' ]
for file_path in file_paths:
   header, endian_test_string, lx1, ly1, nelf, elmap, x, y, vx, vy, vz = read_binary_file(file_path)
   # Append the velocity arrays along with other necessary checks
   vx_list.append(vx)
   vy_list.append(vy)
   vz_list.append(vz)

# Concatenate the velocity arrays along the last dimension (axis=3)
vx = np.concatenate(vx_list, axis=3)
vy = np.concatenate(vy_list, axis=3)
vz = np.concatenate(vz_list, axis=3)

nsteps = vx.shape[3]

# Assuming vx, x, y, nelf are already defined

# Create the figure and axis
fig, ax = plt.subplots()

# Initialize the variables to store vmin and vmax for coloring
vmin = vx[:,:,:,0].min()
vmax = vx[:,:,:,0].max()
dvmin = 1000
dvmax = 0
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

        ax.set_title(f'Plot {j + 1} of {nsteps}')  # Display current plot number (1-based index)

        plt.pause(0.001)