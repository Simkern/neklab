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

def fliplr(f):
    return np.flip(f, axis=1)

def flipudlr(f):
    return np.flip(f, axis=(0,1))

def flipat(f):
    fld = np.empty_like(f)
    a, b = f.shape[:2]
    for i in range(a):
        for j in range(b):
            if len(f.shape)==2:
                fld[i,j]     = f[a-1-j, b-1-i]
            else:
                fld[i,j,...] = f[a-1-j, b-1-i, ...]
    return fld

def fliptr(f):
    if len(f.shape)==2:
        fld = f.T
    else:
        fld = np.empty_like(f)
        for k in range(f.shape[2]):
            fld[:,:,k] = f[:,:,k].T
    return fld

def print_msg(typestr, x1, xt, err):
    print(f'\tx:',' '.join([ f'{d:8.4f}' for d in x1[0,:] ]),' --> ',' '.join([ f'{d:8.4f}' for d in xt[0,:] ]))
    print(f'\t  ',' '.join([ f'{d:8.4f}' for d in x1[1,:] ]),' --> ',' '.join([ f'{d:8.4f}' for d in xt[1,:] ]))
    print(f'error after {typestr}: {err}')

def get_ord(x0,y0,x1,y1,tol=1e-5,debug=False):
    # we assume that the input data contain the flipped data reflecting on which side of the mesh the elements are located
    s0 = x0+y0
    err = abs(s0-(x1+y1)).max()
    if err < tol:
        ord = 0
    if err > tol:
        xtc = fliplr(x1)
        ytc = fliplr(y1)
        err = abs(s0-(xtc+ytc)).max()
        if debug:
            print_msg('lr', x1, xtc, err)
        if err < tol:
            ord = 1
    if err > tol:
        xtc = flipud(x1)
        ytc = flipud(y1)
        err = abs(s0-(xtc+ytc)).max()
        if debug:
            print_msg('ud', x1, xtc, err)
        if err < tol:
            ord = 2
    if err > tol:
        xtc = flipudlr(x1)
        ytc = flipudlr(y1)
        err = abs(s0-(xtc+ytc)).max()
        if debug:
            print_msg('udlr', x1, xtc, err)
        if err < tol:
            ord = 3
    if err > tol:
        xtc = fliptr(x1)
        ytc = fliptr(y1)
        err = abs(s0-(xtc+ytc)).max()
        if debug:
            print_msg('transpose', x1, xtc, err)
        if err < tol:
            ord = 4
    if err > tol:
        xtc = flipat(x1)
        ytc = flipat(y1)
        err = abs(s0-(xtc+ytc)).max()
        if debug:
            print_msg('anti-tranpose', x1, xtc, err)
        if err < tol:
            ord = 5
    if err > tol:
        print(f'The elements could not be mapped!')
        print('Target:')
        print(f'\tx:',' '.join([ f'{d:12.8f}' for d in x0[0,:] ]),' y:',' '.join([ f'{d:12.8f}' for d in y0[0,:] ]))
        print(f'\t  ',' '.join([ f'{d:12.8f}' for d in x0[1,:] ]),'   ',' '.join([ f'{d:12.8f}' for d in y0[1,:] ]))
        print('Origin:')
        print(f'\tx:',' '.join([ f'{d:12.8f}' for d in x1[0,:] ]),' y:',' '.join([ f'{d:12.8f}' for d in y1[0,:] ]))
        print(f'\t  ',' '.join([ f'{d:12.8f}' for d in x1[1,:] ]),'   ',' '.join([ f'{d:12.8f}' for d in y1[1,:] ]))
        print('Attempts:')
        xtc = fliplr(x1); ytc = fliplr(y1)
        err = s0-(xtc+ytc)
        print('LR:')
        print(f'\tx:',' '.join([ f'{d:12.8f}' for d in xtc[0,:] ]),' y:',' '.join([ f'{d:12.8f}' for d in ytc[0,:] ]))
        print(f'\t  ',' '.join([ f'{d:12.8f}' for d in xtc[1,:] ]),'   ',' '.join([ f'{d:12.8f}' for d in ytc[1,:] ]))
        print(f'\terr:',' '.join([ f'{d:12.8e}' for d in err[0,:] ]))
        print(f'\t    ',' '.join([ f'{d:12.8e}' for d in err[1,:] ]))
        print(f'\tError: {abs(err).max():14.12e}')
        xtc = flipud(x1); ytc = flipud(y1)
        err = s0-(xtc+ytc)
        print('UD:')
        print(f'\tx:',' '.join([ f'{d:12.8f}' for d in xtc[0,:] ]),' y:',' '.join([ f'{d:12.8f}' for d in ytc[0,:] ]))
        print(f'\t  ',' '.join([ f'{d:12.8f}' for d in xtc[1,:] ]),'   ',' '.join([ f'{d:12.8f}' for d in ytc[1,:] ]))
        print(f'\terr:',' '.join([ f'{d:12.8e}' for d in err[0,:] ]))
        print(f'\t    ',' '.join([ f'{d:12.8e}' for d in err[1,:] ]))
        print(f'\tError: {abs(err).max():14.12e}')
        xtc = flipudlr(x1); ytc = flipudlr(y1)
        err = s0-(xtc+ytc)
        print('UDLR:')
        print(f'\tx:',' '.join([ f'{d:12.8f}' for d in xtc[0,:] ]),' y:',' '.join([ f'{d:12.8f}' for d in ytc[0,:] ]))
        print(f'\t  ',' '.join([ f'{d:12.8f}' for d in xtc[1,:] ]),'   ',' '.join([ f'{d:12.8f}' for d in ytc[1,:] ]))
        print(f'\terr:',' '.join([ f'{d:12.8e}' for d in err[0,:] ]))
        print(f'\t    ',' '.join([ f'{d:12.8e}' for d in err[1,:] ]))
        print(f'\tError: {abs(err).max():14.12e}')
        xtc = fliptr(x1); ytc = fliptr(y1)
        err = s0-(xtc+ytc)
        print('TR:')
        print(f'\tx:',' '.join([ f'{d:12.8f}' for d in xtc[0,:] ]),' y:',' '.join([ f'{d:12.8f}' for d in ytc[0,:] ]))
        print(f'\t  ',' '.join([ f'{d:12.8f}' for d in xtc[1,:] ]),'   ',' '.join([ f'{d:12.8f}' for d in ytc[1,:] ]))
        print(f'\terr:',' '.join([ f'{d:12.8e}' for d in err[0,:] ]))
        print(f'\t    ',' '.join([ f'{d:12.8e}' for d in err[1,:] ]))
        print(f'\tError: {abs(err).max():14.12e}')
        xtc = flipat(x1); ytc = flipat(y1)
        err = s0-(xtc+ytc)
        print('ATR:')
        print(f'\tx:',' '.join([ f'{d:12.8f}' for d in xtc[0,:] ]),' y:',' '.join([ f'{d:12.8f}' for d in ytc[0,:] ]))
        print(f'\t  ',' '.join([ f'{d:12.8f}' for d in xtc[1,:] ]),'   ',' '.join([ f'{d:12.8f}' for d in ytc[1,:] ]))
        print(f'\terr:',' '.join([ f'{d:12.8e}' for d in err[0,:] ]))
        print(f'\t    ',' '.join([ f'{d:12.8e}' for d in err[1,:] ]))
        print(f'\tError: {abs(err).max():14.12e}')
        sys.exit()
    return ord

def flip_data(f, ord):
    if ord == 0:
        return f
    elif ord == 1: # lr
        return fliplr(f)
    elif ord == 2: # ud
        return flipud(f)
    elif ord == 3: # udlr
        return flipudlr(f)
    elif ord == 4: # tr
        return fliptr(f)
    elif ord == 5: # at
        return flipat(f)
    else:
        print(f'Inconsistent value for ord: {ord}')
        sys.exit()

def get_corners(x,y):
    if x.shape == y.shape:
        if len(x.shape)==2:
            lx1, ly1 = x.shape
            xc = x[[0, 0, lx1-1, lx1-1], [0, ly1-1, 0, ly1-1]].reshape(2,2)
            yc = y[[0, 0, lx1-1, lx1-1], [0, ly1-1, 0, ly1-1]].reshape(2,2)
        elif len(x.shape)==3:
            lx1, ly1, nel = x.shape
            xc = x[[0, 0, lx1-1, lx1-1], [0, ly1-1, 0, ly1-1], :].reshape(2,2,nel)
            yc = y[[0, 0, lx1-1, lx1-1], [0, ly1-1, 0, ly1-1], :].reshape(2,2,nel)
        else:
            print('Input data has the wrong size!')
            print(f'x: {x.shape}')
            print(f'y: {y.shape}')
            sys.exit()
    else:
        print('Input data has unequal size!')
        print(f'x: {x.shape}')
        print(f'y: {y.shape}')
        sys.exit()
    return xc, yc

def symmetrize_fld(dataf, fld, tol=1e-5, debug=False):
    lx1, ly1, nelf = dataf.x.shape

    xs, ys = np.empty_like(dataf.x), np.empty_like(dataf.y)
    fld_s, fld_a = np.empty_like(fld), np.empty_like(fld)

    # get element centers
    xa = np.mean(dataf.x, axis=(0,1))
    ya = np.mean(dataf.y, axis=(0,1))

    # extract corners for the entire mesh
    xc, yc = get_corners(dataf.x, dataf.y)

    # find elements on each side (no cell center can be at 0.0)
    il = [ i for i in range(nelf) if xa[i] < 0.0 ]
    ir = [ i for i in range(nelf) if xa[i] > 0.0 ]
    nl = len(il)
    nr = len(ir)
    if not (nl == nr == int(nelf/2.0)):
        print('Inconsistent mesh split:')
        print(f'Left side:  {nl} elements.')
        print(f'Right side: {nr} elements.')
        print(f'Total: {nelf} elements.')
        sys.exit()

    imap = 0
    mapped = []
    # map elements to left side
    for i in il:
        xl, yl, xlc, ylc, fl = xa[i], ya[i], xc[:,:,i], yc[:,:,i], fld[:,:,i,:]
        found = False
        mindist = 1000.0
        for j in ir:
            xr, yr, xrc, yrc, fr = -xa[j], ya[j], -xc[:,:,j], yc[:,:,j], fld[:,:,j,:]
            d = np.sqrt((xl - xr) ** 2 + (yl - yr) ** 2)
            if d < mindist:
                mindist = d
                imin = j
                x_ = xr
                y_ = yr
            if d < tol:
                order = get_ord(xlc, ylc, xrc, yrc, tol=tol)
                # apply
                fr_l = flip_data(fr, order)
                fs   = (fl + fr_l)/2.0
                fld_s[:,:,i,:] = fs.copy()
                fld_s[:,:,j,:] = flip_data(fs, order)
                fld_a[:,:,i,:] = fl - fld_s[:,:,i,:]
                fld_a[:,:,j,:] = fr - fld_s[:,:,j,:]
                mapped.append([i,j])
                imap += 1
                found = True
        if not found:
            print(f'\nNo match found for element {i} with tolerance {tol:16.12f}')
            print(f'\t centroid: {xl:16.12f}, {yl:16.12f}')
            print(f'Closest match: \n  element {imin} with d = {mindist:16.12f}')
            print(f'\t centroid: {x_:16.12f}, {y_:16.12f}')
            sys.exit()
    if not 2*imap == nelf:
        print(f'imap = {imap}')
        sys.exit()
    return fld_s, fld_a

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

def h2d_to_f2d(h2d_pattern, f2d_file_ref, hfldr='.', outfldr='.', outpattern=None, ord_tol=1e-8):

    datah, metah, nsteps  = read_fields(h2d_pattern, cwd=hfldr)
    dataf, metaf          = read_binary_file(f2d_file_ref, only_mesh=True)

    map_h2f, _ = link_meshes(datah.x, datah.y, dataf.x, dataf.y)

    lx1, ly1, nelfh = datah.vx.shape[:3]
    shp = [lx1, ly1, 2*nelfh, nsteps]
    vx, vy, vz = np.empty(shp), np.empty(shp), np.empty(shp)

    xhc, yhc = get_corners(datah.x, datah.y) # extract corners for the entire mesh
    xfc, yfc = get_corners(dataf.x, dataf.y) # extract corners for the entire mesh

    orientation = [ ]
    sv = [1, -1] # we assume that the first element in each pair in map_h2f is on the left side
    for i, idx in enumerate(map_h2f):
        x0, y0 = xhc[:,:,i], yhc[:,:,i]
        ord = []
        for id, j in enumerate(idx):
            flip = sv[id] # flip the x-axis if the element is on the right side
            x1, y1 = flip*xfc[:,:,j], yfc[:,:,j]
            # get orientation
            order = get_ord(x0, y0, x1, y1, tol=ord_tol)
            # apply
            vx[:,:,j,:] = flip_data(datah.vx[:,:,i,:], order)
            vy[:,:,j,:] = flip_data(datah.vy[:,:,i,:], order)
            vz[:,:,j,:] = flip_data(datah.vz[:,:,i,:], order)
            ord.append(order)
        orientation.append(ord)
    dataf.vx = vx
    dataf.vy = vy
    dataf.vz = vz
    dataf.dt = datah.dt
    metaf.nsave = metah.nsave
    # write file
    if outpattern is None:
        outpattern = h2d_pattern
    
    write_fields(outpattern+'_f', dataf, metaf, nsteps, cwd=outfldr, force=True)

