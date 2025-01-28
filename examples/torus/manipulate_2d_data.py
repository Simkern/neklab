import sys
import numpy as np
import matplotlib.pyplot as plt
import itertools
from matplotlib  import cm
from read_2d_data import read_binary_file, read_fields
from write_2d_data import write_fields
from plot_2d_data import plot_2d_mesh

def flipud(x,y,f):
    return np.flip(x, axis=0), np.flip(y, axis=0), np.flip(f, axis=0)

def antitranspose(x,y,f):
    fld = np.empty_like(f)
    xnew = np.empty_like(x)
    ynew = np.empty_like(y)
    lx1, ly1 = x.shape
    for i in range(lx1):
        for j in range(ly1):
            fld[i,j,:] = f[lx1-1-j, ly1-1-i,:]
            xnew[i,j]  = x[lx1-1-j, ly1-1-i]
            ynew[i,j]  = y[lx1-1-j, ly1-1-i]
    return xnew, ynew, fld

def fliplr(x, y, f):
    return np.flip(x, axis=1), np.flip(y, axis=1), np.flip(f, axis=1)

def symmetrize_fld(x0, y0, fld, if_debug=False):
    lx1, ly1, nelf = x0.shape
    nxy = lx1*ly1
    xa = np.sum(x0, axis=(0,1))/nxy
    ya = np.sum(y0, axis=(0,1))/nxy
    dx = min(np.max(x0, axis=(0,1)) - np.min(x0, axis=(0,1)))
    dy = min(np.max(y0, axis=(0,1)) - np.min(y0, axis=(0,1)))
    tol = (dx+dy)/10e6

    mapped = [ False for i in range(nelf) ]
    imap = 0
    xs, ys = np.empty_like(x0), np.empty_like(x0)
    fld_sym, fld_asym = np.empty_like(fld), np.empty_like(fld)

    for i, (x1, y1) in enumerate(zip(xa, ya)):
        if not mapped[i]:
            for j, (x2, y2) in enumerate(zip(xa, ya)):
                d = np.sqrt((x1 + x2) ** 2 + (y1 - y2) ** 2)
                if d < tol:
                    if x1 > 0:
                        l = i
                        r = j
                    else:
                        l = j
                        r = i
                    idx1 = np.argsort([  x0[a,b,l]+y0[a,b,l] for (a,b) in itertools.product([0,lx1-1], repeat=2) ])
                    idx2 = np.argsort([ -x0[a,b,r]+y0[a,b,r] for (a,b) in itertools.product([0,lx1-1], repeat=2) ])
                    #print(idx2[idx1])
                    s = idx2[idx1]
                    xl   = x0[:,:,l]
                    yl   = y0[:,:,l]
                    fldl = fld[:,:,l,:]
                    if (s == [1, 0, 3, 2]).all():
                        xr, yr, fldr = flipud(x0[:,:,r],y0[:,:,r],fld[:,:,r,:])
                    elif (s == [3, 1, 2, 0]).all():
                        xr, yr, fldr = antitranspose(x0[:,:,r],y0[:,:,r],fld[:,:,r,:])
                    elif (s == [2, 3, 0, 1]).all():
                        xr, yr, fldr = flipud(x0[:,:,r],y0[:,:,r],fld[:,:,r,:])
                    elif (s == [1, 2, 0, 3]).all():
                        xr, yr, fldr = flipud(x0[:,:,r],y0[:,:,r],fld[:,:,r,:])

                    err = abs((xl+yl+xr-yr).max())
                    if err>0.001:
                        if (s == [2, 3, 0, 1]).all():
                            xr, yr, fldr = fliplr(x0[:,:,r],y0[:,:,r],fld[:,:,r,:])
                        if (s == [1, 0, 3, 2]).all():
                            xr, yr, fldr = fliplr(x0[:,:,r],y0[:,:,r],fld[:,:,r,:])
                        elif (s == [1, 2, 0, 3]).all():
                            xr, yr, fldr = fliplr(x0[:,:,r],y0[:,:,r],fld[:,:,r,:])

                    err = abs((xl+yl+xr-yr).max())
                    #print(err)
                    if (err > 0.0001):
                        print(f'x el1:')
                        for el in range(0,lx1,1):
                            print(' '.join([ f' {d:8.4f}' for d in x0[el,range(0,lx1,1),l] ]))
                        print(f'-x el2:')
                        for el in range(0,lx1,1):
                            print(' '.join([ f' {d:8.4f}' for d in -x0[el,range(0,lx1,1),r] ]))
                        print(f'-xr:')
                        for el in range(0,lx1,1):
                            print(' '.join([ f' {d:8.4f}' for d in -xr[el,range(0,lx1,1)] ]))
                        print('\n\n')
                        print('\n\n')
                        sys.exit()
                    flds = (fldl + fldr)/2.0
                    fld_sym[:,:,i] = flds.copy()
                    fld_sym[:,:,j] = flds.copy()
                    fld_asym[:,:,i] = fldl - flds
                    fld_asym[:,:,j] = fldr - flds
                    xs[:,:,i] = xl
                    xs[:,:,j] = xr
                    ys[:,:,i] = yl
                    ys[:,:,j] = yr

                    mapped[j] = True
                    imap += 1
            mapped[i] = True
            imap += 1
    if not imap == nelf:
        print('Not all poitns mapped!')
        sys.exit()
    return xs, ys, fld_sym, fld_asym

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

