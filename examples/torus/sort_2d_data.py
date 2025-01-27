import sys
import numpy as np
from read_2d_data import read_binary_file, read_fields
from write_2d_data import write_fields

def link_meshes(xh, yh, xf, yf, pattern=None):
    
    lx1, ly1, _ = xh.shape
    nxy = lx1*ly1
    xaf = sum(xf, axis=(0,1))/nxy
    yaf = sum(yf, axis=(0,1))/nxy
    xah = sum(xh, axis=(0,1))/nxy
    yah = sum(yh, axis=(0,1))/nxy
    dx = min(np.max(xh, axis=(0,1)) - np.min(xh, axis=(0,1)))
    dy = min(np.max(yh, axis=(0,1)) - np.min(yh, axis=(0,1)))
    
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
            sys.exit()
    
    if map_f2h.min() == 0:
        print('Not all points found!')
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

    write_fields(h2d_pattern, xf, yf, vx, vy, vz, elmapf, dt2d, metadata, nsteps, cwd=ffldr)

