import sys
import numpy as np
import matplotlib.pyplot as plt
import itertools
from matplotlib  import cm
from read_2d_data import read_binary_file, read_fields
from write_2d_data import write_fields
from plot_2d_data import plot_2d_mesh

def flipud(f):
    return np.flip(f, axis=0)

def antitranspose(f):
    fld = np.empty_like(f)
    lx1, ly1 = f.shape[:2]
    for i in range(lx1):
        for j in range(ly1):
            fld[i,j,...] = f[lx1-1-j, ly1-1-i,...]
    return fld

def fliplr(f):
    return np.flip(f, axis=1)

def trans(f):
    return np.flip(f, axis=(0,1))

def symmetrize_fld(x0, y0, fld, if_debug=False):
    lx1, _, nelf = x0.shape
    xa = np.mean(x0, axis=(0,1))
    ya = np.mean(y0, axis=(0,1))
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
                        xr   = flipud(x0[:,:,r])
                        yr   = flipud(y0[:,:,r])
                        fldr = flipud(fld[:,:,r,:])
                    elif (s == [3, 1, 2, 0]).all():
                        xr   = antitranspose(x0[:,:,r])
                        yr   = antitranspose(y0[:,:,r])
                        fldr = antitranspose(fld[:,:,r,:])
                    elif (s == [2, 3, 0, 1]).all():
                        xr   = flipud(x0[:,:,r])
                        yr   = flipud(x0[:,:,r])
                        fldr = flipud(fld[:,:,r,:])
                    elif (s == [1, 2, 0, 3]).all():
                        xr   = flipud(x0[:,:,r])
                        yr   = flipud(y0[:,:,r])
                        fldr = flipud(fld[:,:,r,:])

                    err = abs((xl+yl+xr-yr).max())
                    if err>0.001:
                        if (s == [2, 3, 0, 1]).all():
                            xr   = fliplr(x0[:,:,r])
                            yr   = fliplr(y0[:,:,r])
                            fldr = fliplr(fld[:,:,r,:])
                        if (s == [1, 0, 3, 2]).all():
                            xr   = fliplr(x0[:,:,r])
                            yr   = fliplr(y0[:,:,r])
                            fldr = fliplr(fld[:,:,r,:])
                        elif (s == [1, 2, 0, 3]).all():
                            xr   = fliplr(x0[:,:,r])
                            yr   = fliplr(y0[:,:,r])
                            fldr = fliplr(fld[:,:,r,:])

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
    
    xaf = np.mean(x0f, axis=(0,1))
    yaf = np.mean(y0f, axis=(0,1))
    xah = np.mean(x0h, axis=(0,1))
    yah = np.mean(y0h, axis=(0,1))
    dx = min(np.max(x0h, axis=(0,1)) - np.min(x0h, axis=(0,1)))
    dy = min(np.max(y0h, axis=(0,1)) - np.min(y0h, axis=(0,1)))
    
    tol = (dx+dy)/200

    map_h2f = []
    map_f2h = [ -1 for i in range(len(xaf)) ]
    for i, (xh, yh) in enumerate(zip(xah, yah)):
        idx = []
        xav = []
        for j, (xf, yf) in enumerate(zip(xaf, yaf)):
            d = np.sqrt((xh + np.abs(xf)) ** 2 + (yh - yf) ** 2)
            #print(f'\telement {j:3d}: {xf:8.4f}, {yf:8.4f}, dist = {d:f}')
            if d < tol:
                idx.append(j)
                xav.append(xf)
        if len(idx) == 2:
            # always the left half first
            isort = np.argsort(np.array(xav))
            idx = [ idx[a] for a in isort ]
            map_h2f.append(idx)
            for j in idx:
                map_f2h[j] = i
        else:
            print('Inconsistent meshes!')
            fig, ax = plt.subplots()
            plot_2d_mesh(ax, x0f, y0f, only_edges=True)
            plot_2d_mesh(ax, x0h, y0h, only_edges=True)
            ax.scatter(xaf, yaf, c='k')
            ax.scatter(xah, yah, c='r', marker='x', s=50)
            plt.show()
            sys.exit()
    if -1 in map_f2h:
        print('Not all points found!')
        fig, ax = plt.subplots()
        plot_2d_mesh(ax, x0f, y0f, only_edges=True)
        plot_2d_mesh(ax, x0h, y0h, only_edges=True)
        ax.scatter(xaf, yaf, c='k')
        ax.scatter(xah, yah, c='r', marker='x', s=50)
        idx0 = map_f2h.index(-1)
        print(idx0)
        ax.scatter(xaf[idx0], yaf[idx0], s=200, fc='none', edgecolors='b')
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

    xh, yh, vxh, vyh, vzh, _, dt2d, _, nsteps = read_fields(h2d_pattern, cwd=hfldr)
    xf, yf, _, _, _,       elmapf, _,    metaf = read_binary_file(f2d_file_ref, only_mesh=True)

    map_h2f, _ = link_meshes(xh, yh, xf, yf)

    lx1, ly1, nelfh = vxh.shape[:3]
    shp = [lx1, ly1, 2*nelfh, nsteps]
    vx, vy, vz = np.empty(shp), np.empty(shp), np.empty(shp)

    sv = [1, -1]
    for i, idx in enumerate(map_h2f):
        idx1 = np.argsort([  xh[a,b,i]+yh[a,b,i] for (a,b) in itertools.product([0,lx1-1], repeat=2) ])
        #print(xh[0,:,i])
        #print(f'\nhel {i:3d}: \t{idx1}')
        x0 = xh[:,:,i]
        y0 = yh[:,:,i]
        for id, j in enumerate(idx):
            flip = sv[id]
            idx2 = np.argsort([ flip*xf[a,b,j]+yf[a,b,j] for (a,b) in itertools.product([0,lx1-1], repeat=2) ])
            if (idx1 == idx2).all():
                xt = xf[:,:,j]
                yt = yf[:,:,j]
                vx[:,:,j,:] = vxh[:,:,i,:]
                vy[:,:,j,:] = vyh[:,:,i,:]
                vz[:,:,j,:] = vzh[:,:,i,:]
            else:
                s = idx1[idx2]
                if (s == [0, 1, 2, 3]).all():
                    xt          = xf[:,:,j]
                    yt          = yf[:,:,j]
                    vx[:,:,j,:] = vxh[:,:,i,:]
                    vy[:,:,j,:] = vyh[:,:,i,:]
                    vz[:,:,j,:] = vzh[:,:,i,:]
                elif (s == [3, 2, 1, 0]).all():
                    xt          = trans(xf[:,:,j])
                    yt          = trans(yf[:,:,j])
                    vx[:,:,j,:] = trans(vxh[:,:,i,:])
                    vy[:,:,j,:] = trans(vyh[:,:,i,:])
                    vz[:,:,j,:] = trans(vzh[:,:,i,:])
                elif (s == [3, 1, 2, 0]).all():
                    xt          = antitranspose(xf[:,:,j])
                    yt          = antitranspose(yf[:,:,j])
                    vx[:,:,j,:] = antitranspose(vxh[:,:,i,:])
                    vy[:,:,j,:] = antitranspose(vyh[:,:,i,:])
                    vz[:,:,j,:] = antitranspose(vzh[:,:,i,:])
                elif (s == [1, 0, 3, 2]).all():
                    xt          = fliplr(xf[:,:,j])
                    yt          = fliplr(yf[:,:,j])
                    vx[:,:,j,:] = fliplr(vxh[:,:,i,:])
                    vy[:,:,j,:] = fliplr(vyh[:,:,i,:])
                    vz[:,:,j,:] = fliplr(vzh[:,:,i,:])
                elif (s == [0, 2, 1, 3]).all():
                    xt              = xf[:,:,j].T
                    yt              = yf[:,:,j].T
                    for l in range(vx.shape[-1]):
                        vx[:,:,j,l] = vxh[:,:,i,l].T
                        vy[:,:,j,l] = vyh[:,:,i,l].T
                        vz[:,:,j,l] = vzh[:,:,i,l].T
            err = abs((x0+y0-flip*xt-yt).max())
            if (err > 0.0001):
                if (s == [3, 1, 2, 0]).all():
                    xt              = xf[:,:,j].T
                    yt              = yf[:,:,j].T
                    for l in range(vx.shape[-1]):
                        vx[:,:,j,l] = vxh[:,:,i,l].T
                        vy[:,:,j,l] = vyh[:,:,i,l].T
                        vz[:,:,j,l] = vzh[:,:,i,l].T
                elif (s == [0, 3, 2, 1]).all():
                    xt              = xf[:,:,j].T
                    yt              = yf[:,:,j].T
                    for l in range(vx.shape[-1]):
                        vx[:,:,j,l] = vxh[:,:,i,l].T
                        vy[:,:,j,l] = vyh[:,:,i,l].T
                        vz[:,:,j,l] = vzh[:,:,i,l].T
                elif (s == [0, 1, 2, 3]).all():
                    xt          = trans(xf[:,:,j])
                    yt          = trans(yf[:,:,j])
                    vx[:,:,j,:] = trans(vxh[:,:,i,:])
                    vy[:,:,j,:] = trans(vyh[:,:,i,:])
                    vz[:,:,j,:] = trans(vzh[:,:,i,:])
                elif (s == [2, 1, 3, 0]).all():
                    xt          = trans(xf[:,:,j])
                    yt          = trans(yf[:,:,j])
                    vx[:,:,j,:] = trans(vxh[:,:,i,:])
                    vy[:,:,j,:] = trans(vyh[:,:,i,:])
                    vz[:,:,j,:] = trans(vzh[:,:,i,:])
                elif (s == [0, 2, 1, 3]).all():
                    xt          = antitranspose(xf[:,:,j])
                    yt          = antitranspose(yf[:,:,j])
                    vx[:,:,j,:] = antitranspose(vxh[:,:,i,:])
                    vy[:,:,j,:] = antitranspose(vyh[:,:,i,:])
                    vz[:,:,j,:] = antitranspose(vzh[:,:,i,:])
                err = abs((x0+y0-flip*xt-yt).max())
                if (err > 0.0001):
                    print(f'Error = {err}')
                    print(f'       : \t', end='')
                    print(' '.join([ f'{d:8.4f}' for d in [ xh[x,y,i] for (x,y) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1])] ]))
                    print(f'       : \t', end='')
                    print(' '.join([ f'{d:8.4f}' for d in [ yh[x,y,i] for (x,y) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1])] ]))
                    print(s)
                    print(f'\tfel {j:3d}:')
                    print(f'\t       :', end='')
                    print(' '.join([ f'{d:8.4f}' for d in [ flip*xf[x,y,j] for (x,y) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1])] ]))
                    print(f'\t       :', end='')
                    print(' '.join([ f'{d:8.4f}' for d in [ yf[x,y,j] for (x,y) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1])] ]))
                    print(f'\tflip{j:3d}:', end='')
                    print(' '.join([ f'{d:8.4f}' for d in [ flip*xt[x,y] for (x,y) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1])] ]))
                    print(f'\tflip{j:3d}:', end='')
                    print(' '.join([ f'{d:8.4f}' for d in [ yt[x,y] for (x,y) in zip([0,0,lx1-1,lx1-1],[0,lx1-1,0,lx1-1])] ]))
                    sys.exit()
    write_fields(h2d_pattern+'_f', xf, yf, vx, vy, vz, elmapf, dt2d, metaf, nsteps, cwd=ffldr, force=True)

