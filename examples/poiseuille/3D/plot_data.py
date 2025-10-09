import sys, os
import numpy as np
import pymech as pm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def read_his(filename, nslice=1):
    print(f'Number of slices: {nslice}')
    with open(filename, 'r') as f:
        # Read number of points
        n = int(f.readline().strip())
        print(f'Number of points: {n}')

        # Read next n lines: coordinates
        coords_raw = np.loadtxt([next(f) for _ in range(n)])
        print(f'Read coordinates.')

        # Read next n lines: time, ux, uy, uz, p
        data_raw = np.loadtxt(f)
        print(f'Read data.')

    # Reshape data: (m, nslice, n, 5)
    if not np.mod(n, nslice) == 0:
        print('Number of points not divisible by input number of slices.')
    npts_slice = n // nslice
    nstep = data_raw.shape[0] // npts_slice // nslice
    data = data_raw.reshape((nstep, nslice, npts_slice, 5))
    coords = coords_raw.reshape((nslice, npts_slice, 3))
    print(f'Coords: {coords.shape}')
    print(f'Data:   {data.shape}')
    z_list = [ coords[i,0,0] for i in range(nslice) ]

    return coords, data, z_list

def write_data(coords, data, output_dir="."):
    """
    Writes coordinates and data to separate files.
    
    Parameters:
    - coords: numpy array of shape (n, 3)
    - data: numpy array of shape (m, n, 5)
    - output_dir: directory where files will be saved
    """
    # Write coordinates
    coord_path = f"{output_dir}/coords.dat"
    np.savetxt(coord_path, coords, fmt="%.15e")
    print(f"Wrote coordinates to {coord_path}")

    # Write each time instance
    for t in range(data.shape[0]):
        file_path = f"{output_dir}/data_{t:03d}.dat"
        np.savetxt(file_path, data[t], fmt="%.15e")
        print(f"Wrote time step {t} data to {file_path}")

def load_coords_and_data(istep=-1):
    coords = np.loadtxt("coords.dat")
    data_files = sorted([
        f for f in os.listdir('.')
        if f.startswith("data_") and f.endswith(".dat")
    ])
    if istep >= 0:
        data = [np.loadtxt(f) for f in [data_files[istep]]]
    else:
        data = [np.loadtxt(f) for f in data_files]
    return coords, data

def generate_mesh(coords):
    xy, idx = np.unique(np.stack((coords[0, :, 0], coords[0, :, 1])).T, axis=0, return_index=True)
    xmin, xmax = xy[:,0].min(), xy[:,0].max()
    ymin, ymax = xy[:,1].min(), xy[:,1].max()
    xx, yy = np.meshgrid(np.linspace(xmin, xmax, num=200, endpoint=True), np.linspace(ymin, ymax, num=200, endpoint=True))
    xy_plot = np.stack([xx.ravel(), yy.ravel()], axis=1)
    return xy, xy_plot, xx, yy, idx

def plot_t(fig, ax, xy, xy_plot, xx, yy, idx, data_list, nfield_index=4, istep=0, islice=0):
    fld_names = [ 'u_x', 'u_y', 'u_z', 'pressure' ]
    # get data
    field = data_list[istep][islice, :, nfield_index]
    # interpolate
    fld = griddata(xy, field[idx], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
    # plot
    if islice == 0:
        comp = 'Re'
    else:
        comp = 'Im'
    c = ax.contourf(xx, yy, fld, levels=50, cmap='viridis')
    ax.set_title(f"Instant {istep}, fld: {fld_names[field_index]}, {comp}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(c, ax=ax, orientation='vertical')

def plot_all(xy, xy_plot, xx, yy, idx, data_list, field_index=4, islice=0):
    fld_names = [ 'u_x', 'u_y', 'u_z', 'pressure' ]
    fig, axs = plt.subplots(5, 3, figsize=(10, 15), constrained_layout=True)
    axs = axs.flatten()

    for istep, data in enumerate(data_list):
        # get data
        field = data[islice, :, field_index]
        # interpolate
        fld = griddata(xy, field[idx], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
        # plot
        ax = axs[istep]
        if islice == 0:
           comp = 'Re'
        else:
           comp = 'Im'
        c = ax.contourf(xx, yy, fld, levels=50, cmap='viridis')
        ax.set_title(f"Instant {istep}, fld: {fld_names[field_index]}, {comp}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        # Optional: add colorbar
        fig.colorbar(c, ax=ax, orientation='vertical')

    # Hide any unused subplots (shouldn't happen if 15 time steps)
    for j in range(len(data_list), len(axs)):
        axs[j].axis('off')
    plt.show()

# Write them to disk
#coords, data = read_his("poiseuille.his")
#print("Coordinates shape:", coords.shape)  # (n, 3)
#print("Data shape:", data.shape)          # (15, n, 5)
#write_data(coords, data)

nslice = 2
coords, data_list, z_list = read_his("poiseuille.his", nslice=nslice)
xy, xy_plot, xx, yy, idx = generate_mesh(coords)

basename = 'poiseuille0.f00001'
fldr2Dh = '../2Dh/'

#2Dh nek data
fld_names = [ 'u_x', 'u_y', 'u_z', 'pressure' ]
fig, axs = plt.subplots(4,4, figsize=(10, 15), constrained_layout=True)
for jp in range(2):
    fname = os.path.join(fldr2Dh, 'rv'+str(jp+1)+basename)
    print(f'Read {fname}')
    m2d = pm.neksuite.readnek(fname)
    x2d = np.array([ e.pos[0].ravel() for e in m2d.elem ]).ravel()
    y2d = np.array([ e.pos[1].ravel() for e in m2d.elem ]).ravel()
    xy_nek = np.stack((x2d[idx], y2d[idx])).T

    if jp == 0:
       comp = 'Re'
    else:
       comp = 'Im'

    u = np.array([ e.vel[0].ravel() for e in m2d.elem ]).ravel()
    v = np.array([ e.vel[1].ravel() for e in m2d.elem ]).ravel()
    w = np.array([ e.temp[0].ravel() for e in m2d.elem ]).ravel()
    p = np.array([ e.pres[0].ravel() for e in m2d.elem ]).ravel()

    flds = [u, v, w, p]

    iplt = jp
    for ifld in range(4):
      ax = axs[ifld,iplt]
      fld = griddata(xy_nek, flds[ifld][idx], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
      c = ax.contourf(xx, yy, fld, levels=50, cmap='viridis')
      ax.set_title(f"Instant {0}, fld: {fld_names[ifld]}, {comp}")
      fig.colorbar(c, ax=ax, orientation='vertical')
      ax.set_xlabel("x")
      ax.set_ylabel("y")
    #fig.suptitle('residual velocity')

for jp in range(2):
    fname = os.path.join(fldr2Dh, 'dv'+str(jp+1)+basename)
    print(f'Read {fname}')
    m2d = pm.neksuite.readnek(fname)
    x2d = np.array([ e.pos[0].ravel() for e in m2d.elem ]).ravel()
    y2d = np.array([ e.pos[1].ravel() for e in m2d.elem ]).ravel()
    xy_nek = np.stack((x2d[idx], y2d[idx])).T

    if jp == 0:
       comp = 'Re'
    else:
       comp = 'Im'

    u = np.array([ e.vel[0].ravel() for e in m2d.elem ]).ravel()
    v = np.array([ e.vel[1].ravel() for e in m2d.elem ]).ravel()
    w = np.array([ e.temp[0].ravel() for e in m2d.elem ]).ravel()
    p = np.array([ e.pres[0].ravel() for e in m2d.elem ]).ravel()

    flds = [u, v, w, p]

    iplt = 2 + jp
    for ifld in range(4):
      ax = axs[ifld,iplt]
      fld = griddata(xy_nek, flds[ifld][idx], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
      c = ax.contourf(xx, yy, fld, levels=50, cmap='viridis')
      ax.set_title(f"Instant {1}, fld: {fld_names[ifld]}, {comp}")
      fig.colorbar(c, ax=ax, orientation='vertical')
      ax.set_xlabel("x")
      ax.set_ylabel("y")
fig.suptitle('velocity 2Dh')









fld_names = [ 'pressure' ]
fig, axs = plt.subplots(1,4, figsize=(10, 4), constrained_layout=True)
for jp in range(2):
    fname = os.path.join(fldr2Dh, 'dp'+str(jp+1)+basename)
    print(f'Read {fname}')
    m2d = pm.neksuite.readnek(fname)
    x2d = np.array([ e.pos[0].ravel() for e in m2d.elem ]).ravel()
    y2d = np.array([ e.pos[1].ravel() for e in m2d.elem ]).ravel()
    xy_nek = np.stack((x2d[idx], y2d[idx])).T

    if jp == 0:
       comp = 'Re'
    else:
       comp = 'Im'

    p = np.array([ e.pres[0].ravel() for e in m2d.elem ]).ravel()

    flds = [p]

    iplt = jp
    ifld = 0
    ax = axs[iplt]
    fld = griddata(xy_nek, flds[ifld][idx], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
    c = ax.contourf(xx, yy, fld, levels=50, cmap='viridis')
    ax.set_title(f"Instant {2}, fld: {fld_names[ifld]}, {comp}")
    fig.colorbar(c, ax=ax, orientation='vertical')
    ax.set_xlabel("x")
    ax.set_ylabel("y")

for jp in range(2):
    fname = os.path.join(fldr2Dh, 'ps'+str(jp+1)+basename)
    print(f'Read {fname}')
    m2d = pm.neksuite.readnek(fname)
    x2d = np.array([ e.pos[0].ravel() for e in m2d.elem ]).ravel()
    y2d = np.array([ e.pos[1].ravel() for e in m2d.elem ]).ravel()
    xy_nek = np.stack((x2d[idx], y2d[idx])).T

    if jp == 0:
       comp = 'Re'
    else:
       comp = 'Im'

    p = np.array([ e.pres[0].ravel() for e in m2d.elem ]).ravel()

    flds = [p]

    iplt = 2 + jp
    ifld = 0
    ax = axs[iplt]
    fld = griddata(xy_nek, flds[ifld][idx], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
    c = ax.contourf(xx, yy, fld, levels=50, cmap='viridis')
    ax.set_title(f"Instant {3}, fld: {fld_names[ifld]}, {comp}")
    fig.colorbar(c, ax=ax, orientation='vertical')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
fig.suptitle('pressure 2Dh')




print('\nlogfile:')
with open('logfile.txt', 'r') as f:
    for line in f.readlines():
        if 'HPTS' in line:
            print(line, end='')
print('')



fig, axs = plt.subplots(4,4, figsize=(10, 15), constrained_layout=True)
for islice in range(nslice):
   for field_index in range(4):
       plot_t(fig, axs[field_index,   islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=0, islice=islice)
       plot_t(fig, axs[field_index, 2+islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=1, islice=islice)
fig.suptitle('velocity 3D')

fig, axs = plt.subplots(1,4, figsize=(10, 4), constrained_layout=True)
for islice in range(nslice):
   field_index = 3
   plot_t(fig, axs[  islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=2, islice=islice)
   plot_t(fig, axs[2+islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=3, islice=islice)
fig.suptitle('pressure 3D')


plt.show()

