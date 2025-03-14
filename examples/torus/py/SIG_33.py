import numpy as np
import matplotlib.pyplot as plt

from read_2d_data import read_fields
from plot_2d_data import plot_2d_fld

pattern = 'prod/run/Wo_035.0/Q_0.280/n2dtorus'

data, meta, nsteps = read_fields(pattern, only_mesh=False)

# Create the figure and axis
fig, axs = plt.subplots(2, 3, figsize=(15, 15))

nplots = np.prod(axs.shape)
T = np.sum(data.dt)

istp = [ int(nsteps*i/nplots) for i in np.arange(nplots+1) if i < nplots ]

for ip, ax in enumerate(axs.flatten()):
    istep = istp[ip]
    plot_2d_fld(ax, data.x, data.y, data.vx, istep=istep, draw_elements=False, draw_mesh=False, if_half=True)
    t = np.sum(data.dt[:istep])/T*100
    ax.set_title(f'{t:3.0f}'+r' $\% \:T$', fontsize=30)  # Display current plot number (1-based index)
    ax.axis('off')
plt.show()