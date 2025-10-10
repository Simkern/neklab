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

    # Reshape data: (nstep, nslice, n, 5)
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

def extract_nek(fname, names, jp):
    print(f'Read {fname}')
    m2d = pm.neksuite.readnek(fname)
    x2d = np.array([ e.pos[0].ravel() for e in m2d.elem ]).ravel()
    y2d = np.array([ e.pos[1].ravel() for e in m2d.elem ]).ravel()
    xy_nek = np.stack((x2d[idx], y2d[idx])).T

    if jp == 0:
       comp = 'Re'
    else:
       comp = 'Im'
 
    flds = []
    if 'u_x' in names:
        fld = np.array([ e.vel[0].ravel() for e in m2d.elem ]).ravel()[idx]
        print(f'\tu_x: min/max:  {fld.min():.8f}/{fld.min():.8f}')
        flds.append(fld)
    if 'u_y' in names:
        fld = np.array([ e.vel[1].ravel() for e in m2d.elem ]).ravel()[idx]
        print(f'\tu_y: min/max:  {fld.min():.8f}/{fld.min():.8f}')
        flds.append(fld)
    if 'u_z' in names:
        fld = np.array([ e.temp[0].ravel() for e in m2d.elem ]).ravel()[idx]
        print(f'\tu_z: min/max:  {fld.min():.8f}/{fld.min():.8f}')
        flds.append(fld)
    if 'pr' in names:
        fld = np.array([ e.pres[0].ravel() for e in m2d.elem ]).ravel()[idx]
        print(f'\tpr:  min/max:  {fld.min():.8f}/{fld.min():.8f}')
        flds.append(fld)

    return flds, xy_nek, comp

def plot_t(fig, ax, xy, xy_plot, xx, yy, idx, data_list, field_index=4, istep=0, islice=0, names=None):
    fld_names = [ 'u_x', 'u_y', 'u_z', 'pr' ]
    name = fld_names[field_index-1]
    # get data
    field = data_list[istep][islice, :, field_index]
    # interpolate
    fld = griddata(xy, field[idx], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
    print(f'\tistep {istep:2d} 3D: {name:8s}: min/max:  {fld.min():.8f}/{fld.min():.8f}')
    # plot
    if islice == 0:
        comp = 'Re'
    else:
        comp = 'Im'
    c = ax.contourf(xx, yy, fld, levels=50, cmap='viridis')
    tle = f'Instant {istep}'
    if names is not None:
        tle = names[istep]
    #ax.set_title(f"{tle}: {name}, {comp}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(c, ax=ax, orientation='vertical')

def plot_all(xy, xy_plot, xx, yy, idx, data_list, field_index=4, islice=0, names=None):
    fld_names = [ 'u_x', 'u_y', 'u_z', 'pressure' ]
    name = fld_names[field_index-1]
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
        #ax.set_title(f"Instant {istep}, fld: {name}, {comp}")
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
plt.close('all')

ifplot = True

nslice = 2
coords, data_list, z_list = read_his("poiseuille.his", nslice=nslice)
xy, xy_plot, xx, yy, idx = generate_mesh(coords)

print('\nlogfile:')
istep = 0
with open('logfile.txt', 'r') as f:
    for line in f.readlines():
        if 'HPTS' in line:
            print(f'{istep:3d}: {line.strip()}')
            istep += 1
print('')
hpts_names_ref = ['resv','dv','resp','dp','dvdp','v','p']
fld_names_ref = [ 'u_x', 'u_y', 'u_z', 'pr' ]

basename = 'poiseuille0.f'
fldr2Dh = '../2Dh/'

#2Dh nek data
steps     = [ 1,2 ]
hpts_names = [name + str(step) for step in steps for name in hpts_names_ref ]
fld_nms   = [ [ 'u_x', 'u_y', 'u_z' ], 
              [ 'pr' ], 
              [ 'u_x', 'u_y', 'u_z' ],
              [ 'u_x', 'u_y', 'u_z' ], 
              [ 'pr' ] ]
fld_tle   = [ [ 'resv', 'dv' ],
              [ 'resp', 'dp' ],
              [ 'dvdp' ], 
              [ 'v' ],
              [ 'p' ] ]
fld_pfx   = [ [ 'rv', 'dv' ], 
              [ 'rp', 'dp' ], 
              [ 'vv' ],
              [ 'vl' ],
              [ 'pr' ] ]

#fld_nms   = [ [ 'u_z' ] ]
#fld_tle   = [ [ 'resv', 'dv' ] ] 
#fld_tle   = [ [ 'dvdp' ] ]
#fld_pfx   = [ [ 'rv', 'dv' ] ]
#fld_pfx   = [ [ 'vv' ] ]

nplot = sum(len(a) * len(b) for a, b in zip(fld_pfx, fld_nms))

for istep in steps:  # timesteps
    #print(f'Timestep {istep}')
    
    if ifplot:
        fig, axs = plt.subplots(6,nplot, figsize=(18, 8), constrained_layout=True)

    for jp in range(2):
        #print(f'  jp = {jp}')

        icol = 0
        for (fld_names, pfx_names, titles) in zip(fld_nms, fld_pfx, fld_tle):
            fld_idxs  = [ fld_names_ref.index(iname)+1 for iname in fld_names ]  # +1 to skip time

            irow = 3*jp
            for (fld_idx, name) in zip(fld_idxs, fld_names):
                print(f'    {name}')
                for icase, (pfx, tle) in enumerate(zip(pfx_names, titles)):
                    fname = os.path.join(fldr2Dh, pfx+str(jp+1)+basename+f'{istep:05d}')
                    flds, xy_nek, comp = extract_nek(fname, [name], jp)
                    #print(f'      {pfx} {tle}: slice {jp}')

                    if ifplot:
                        ax = axs[irow+0,icol]
                        fld_2D = griddata(xy_nek, flds[0], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
                        c = ax.contourf(xx, yy, fld_2D, levels=50, cmap='viridis')
                        fig.colorbar(c, ax=ax, orientation='vertical')
                        if jp == 0:
                            ax.set_title(f'{name} {tle}', fontsize=15)

                        i3D = hpts_names.index(tle+str(istep)) # get index in hpts list ordering 3dstep = 1 outpost
                        ax = axs[irow+1, icol]
                        field = data_list[i3D][jp, :, fld_idx]
                        # interpolate
                        fld_3D = griddata(xy, field[idx], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
                        c = ax.contourf(xx, yy, fld_3D, levels=50, cmap='viridis')
                        fig.colorbar(c, ax=ax, orientation='vertical')
                        ax = axs[irow+2, icol]
                        c = ax.contourf(xx, yy, fld_2D-fld_3D, levels=50, cmap='RdBu_r')
                        fig.colorbar(c, ax=ax, orientation='vertical')
                        icol += 1
                        #plot_t(fig, ax, xy, xy_plot, xx, yy, idx, data_list, field_index=fld_idx, istep=i3D, islice=jp)
    if ifplot:
        for ax in axs.flatten():
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        fig.suptitle(f'step {istep}', fontsize=25)
        axs[0,0].text(-1, 0, rf'2Dh  real', fontsize=15, rotation=90, va='center', ha='center')
        axs[1,0].text(-1, 0, rf'3D   real', fontsize=15, rotation=90, va='center', ha='center')
        axs[2,0].text(-1, 0, rf'diff real', fontsize=15, rotation=90, va='center', ha='center')
        axs[3,0].text(-1, 0, rf'2Dh  imag', fontsize=15, rotation=90, va='center', ha='center')
        axs[4,0].text(-1, 0, rf'3D   imag', fontsize=15, rotation=90, va='center', ha='center')
        axs[5,0].text(-1, 0, rf'duff imag', fontsize=15, rotation=90, va='center', ha='center')


















'''
for istep in steps:
    if ifplot:
        fig, axs = plt.subplots(len(fld_names),4, figsize=(10, 15), constrained_layout=True)
    icol = 0
    for icase, (pfx, tle) in enumerate(zip(pfx_names, titles)):
        for jp in range(2):
            fname = os.path.join(fldr2Dh, pfx+str(jp+1)+basename+f'{istep:05d}')
            flds, xy_nek, comp = extract_nek(fname, fld_names, jp)

            if ifplot:
                iplt = icol + jp
                for ifld, (fld, name) in enumerate(zip(flds, fld_names)):
                    ax = axs[ifld,iplt]
                    fld = griddata(xy_nek, fld, xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
                    c = ax.contourf(xx, yy, fld, levels=50, cmap='viridis')
                    ax.set_title(f"{tle:6s}{ifld+1}, fld: {name}, {comp}")
                    fig.colorbar(c, ax=ax, orientation='vertical')
        icol += 2

    if ifplot:
        for ax in axs.flatten():
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        fig.suptitle(f'velocity 2Dh, step {istep}')
'''

'''
fig, axs = plt.subplots(3,4, figsize=(10, 15), constrained_layout=True)
for islice in range(nslice):
    for field_index in range(3):
        plot_t(fig, axs[field_index,   islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=0, islice=islice)
        plot_t(fig, axs[field_index, 2+islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=1, islice=islice)
fig.suptitle('velocity 3D')



fld_names = [ 'pr' ]
pfx_names = [ 'r1', 'rp', 'dp' ]
titles    = [ 'div2D', 'div3D', 'dp' ]
steps     = [ 1, 2 ]
poffset = 0
for istep in steps:
    fig, axs = plt.subplots(1, len(pfx_names)*2, figsize=(15, 4), constrained_layout=True)
    iplt = 0
    for ipfx, (pfx, tle) in enumerate(zip(pfx_names,titles)):
        for jp in range(2):
            fname = os.path.join(fldr2Dh, pfx+str(jp+1)+basename+f'{istep:05d}')
            flds, xy_nek, comp = extract_nek(fname, fld_names, jp)

            ax = axs[iplt]
            fld = griddata(xy_nek, flds[0], xy_plot, method='cubic', fill_value=0).reshape(xx.shape)
            c = ax.contourf(xx, yy, fld, levels=50, cmap='viridis')
            ax.set_title(f"{tle:6s}, fld: {fld_names[0]}, {comp}")
            fig.colorbar(c, ax=ax, orientation='vertical')
            iplt += 1

    for ax in axs.flatten():
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.suptitle(f'pr 2Dh, step {istep}')








fig, axs = plt.subplots(1,4, figsize=(10, 4), constrained_layout=True)
for islice in range(nslice):
    field_index = 3
    plot_t(fig, axs[  islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=2, islice=islice, names=hpts_names)
    plot_t(fig, axs[2+islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=3, islice=islice, names=hpts_names)
fig.suptitle('pressure 3D, step 1')

fig, axs = plt.subplots(1,4, figsize=(10, 4), constrained_layout=True)
for islice in range(nslice):
    field_index = 3
    plot_t(fig, axs[  islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=6, islice=islice, names=hpts_names)
    plot_t(fig, axs[2+islice], xy, xy_plot, xx, yy, idx, data_list, nfield_index=field_index+1, istep=7, islice=islice, names=hpts_names)
fig.suptitle('pressure 3D, step 2')
'''


plt.show()