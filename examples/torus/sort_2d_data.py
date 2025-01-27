import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib  import cm
from read_2d_data import read_binary_file, read_fields
from write_2d_data import write_fields
from plot_2d_data import plot_2d_mesh

def symmetrize_fld(x0, y0, fld):
    lx1, ly1, npts = x0.shape
    nxy = lx1*ly1
    xa = np.sum(x0, axis=(0,1))/nxy
    ya = np.sum(y0, axis=(0,1))/nxy
    dx = min(np.max(x0, axis=(0,1)) - np.min(x0, axis=(0,1)))
    dy = min(np.max(y0, axis=(0,1)) - np.min(y0, axis=(0,1)))
    tol = (dx+dy)/200

    unique = []
    mapped = [ False for i in range(npts) ]

    for i, (x1, y1) in enumerate(zip(xa, ya)):
        if not mapped[i]:
            for j, (x2, y2) in enumerate(zip(xa, ya)):
                d = np.sqrt((x1 + x2) ** 2 + (y1 - y2) ** 2)
                if d < tol:
                    if x1 > 0:
                        left = i
                        right = j
                    else:
                        left = j
                        right = i
                    print(f'                ', end='')
                    for (r,c) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1]):
                        print(f'   {r} {c}   ', end='')
                    print(f'')
                    print(f'element {i+1:3d}: rx ', end='')
                    for (r,c) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1]):
                        print(f' {x0[r,c,right]:8.5f}', end='')
                    print(f'')
                    print(f'           : ry ', end='')
                    for (r,c) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1]):
                        print(f' {y0[r,c,right]:8.5f}', end='')
                    print(f'')
                    print(f'           : lx ', end='')
                    for (r,c) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1]):
                        print(f' {x0[r,c,left]:8.5f}', end='')
                    print(f'')
                    print(f'           : ly ', end='')
                    for (r,c) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1]):
                        print(f' {y0[r,c,left]:8.5f}', end='')
                    print(f'')
                    print(f'                ', end='')
                    for (r,c) in zip([lx1-1,lx1-1,0,0],[0,lx1-1,0,lx1-1]):
                        print(f'   {r} {c}   ', end='')
                    print(f'')
                    print(f'flipped    : lx ', end='')
                    for (r,c) in zip([lx1-1,lx1-1,0,0],[0,lx1-1,0,lx1-1]):
                        print(f' {x0[r,c,left]:8.5f}', end='')
                    print(f'')
                    print(f'           : ly ', end='')
                    for (r,c) in zip([lx1-1,lx1-1,0,0],[0,lx1-1,0,lx1-1]):
                        print(f' {y0[r,c,left]:8.5f}', end='')
                    print(f'\n')
                    unique.append([i,j])
    
    fld_sym, fld_asym = np.empty_like(fld), np.empty_like(fld)
    for i, (el1, el2) in enumerate(unique):
        fld_sym[:,:,i,:] = (fld[:,:,el1,:] + np.flip(fld[:,:,el2,:], axis=(0,1)))/2.0
        fld_asym[:,:,i,:] = fld[:,:,i,:] - fld_sym[:,:,i,:]
    
    return fld_sym, fld_asym

def link_meshes(x0h, y0h, x0f, y0f, pattern=None):
    
    lx1, ly1, _ = x0h.shape
    nxy = lx1*ly1
    xaf = np.sum(x0f, axis=(0,1))/nxy
    yaf = np.sum(y0f, axis=(0,1))/nxy
    xah = np.sum(x0h, axis=(0,1))/nxy
    yah = np.sum(y0h, axis=(0,1))/nxy
    dx = min(np.max(x0h, axis=(0,1)) - np.min(x0h, axis=(0,1)))
    dy = min(np.max(y0h, axis=(0,1)) - np.min(y0h, axis=(0,1)))
    
    tol = (dx+dy)/200

    map_h2f = []
    map_f2h = [ 0 for i in range(len(xaf)) ]
    for i, (xh, yh) in enumerate(zip(xah, yah)):
        idx = []
        for j, (xf, yf) in enumerate(zip(xaf, yaf)):
            d = np.sqrt((xh - xf) ** 2 + (yh - yf) ** 2)
            if d < tol:
                idx.append(j)
        if len(idx) == 2:
            map_h2f.append(idx)
            map_f2h[idx] = i
        else:
            print('Inconsistent meshes!')
            fig, ax = plt.subplots()
            plot_2d_mesh(ax, x0f, y0f, only_edges=True)
            plot_2d_mesh(ax, x0h, y0h, only_edges=True)
            ax.scatter(xaf, yaf, c='k')
            ax.scatter(xah, yah, c='r', marker='x', s=50)
            plt.show()
            sys.exit()
    
    if map_f2h.min() == 0:
        print('Not all points found!')
        fig, ax = plt.subplots()
        plot_2d_mesh(ax, x0f, y0f, only_edges=True)
        plot_2d_mesh(ax, x0h, y0h, only_edges=True)
        ax.scatter(xaf, yaf, c='k')
        ax.scatter(xah, yah, c='r', marker='x', s=50)
        plt.show()
        sys.exit()

    if pattern is not None:
        with open(pattern+'_h2f.txt', 'w') as f:
            for i, idx in enumerate(map_h2f):
                f.write(f'{i}: {idx}')
        with open(pattern+'_f2h.txt', 'w') as f:
            for i, j in enumerate(map_f2h):
                f.write(f'{i}: {j}')

    return map_h2f, map_f2h

def h2d_to_f2d(h2d_pattern, hfldr, f2d_file_ref, ffldr):

    xh, yh, vxh, vyh, vzh, _, dt2d, metadata, nsteps = read_fields(h2d_pattern, cwd=hfldr)
    xf, yf, _, _, _, elmapf, _, _ = read_binary_file(f2d_file_ref, only_mesh=True)

    _, map_f2h = link_meshes(xh, yh, xf, yf)

    lx1, ly1, nelfh = vxh.shape[:3]
    vx, vy, vz = np.empty([lx1, ly1, 2*nelfh, nsteps])

    for i, j in enumerate(map_f2h):
        vx[:,:,i,:] = vxh[:,:,j,:]
        vy[:,:,i,:] = vyh[:,:,j,:]
        vz[:,:,i,:] = vzh[:,:,j,:]

    write_fields(h2d_pattern+'_f', xf, yf, vx, vy, vz, elmapf, dt2d, metadata, nsteps, cwd=ffldr)

