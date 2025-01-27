import sys, os
import struct
import glob
import numpy as np

def read_int(f,emode,nvar):
    """read integer array"""
    isize = 4
    llist = f.read(isize*nvar)
    llist = list(struct.unpack(emode+nvar*'i', llist))
    return llist

def read_flt(f,emode,wdsize,nvar):
    """read real array"""
    if (wdsize == 4):
        realtype = 'f'
    elif (wdsize == 8):
        realtype = 'd'
    llist = f.read(wdsize*nvar)
    llist = np.frombuffer(llist, dtype=emode+realtype, count=nvar)
    return llist

def read_fields(filepattern, only_mesh=False, cwd='.'):

    pat = os.path.join(cwd, filepattern+'*')

    files = sorted(glob.glob(pat))
    # Check if there are any files that match the pattern
    nfiles = len(files)
    if files:
        print(f"Found {nfiles} file(s):")
        for file in files:
            print(f'   {file}')
    else:
        print(f"No files found matching the pattern '{filepattern}'.")
        sys.exit()

    vx_list, vy_list, vz_list, dt2d_list = [], [], [], []
    for file in files:
        x, y, vxr, vyr, vzr, elmap, dt2d, metadata = read_binary_file(file, only_mesh)
        # Append the velocity arrays along with other necessary checks
        vx_list.append(vxr)
        vy_list.append(vyr)
        vz_list.append(vzr)
        dt2d_list.append(dt2d)

    # Concatenate the velocity arrays along the last dimension (axis=3)
    vx = np.concatenate(vx_list, axis=3)
    vy = np.concatenate(vy_list, axis=3)
    vz = np.concatenate(vz_list, axis=3)
    dt2d = np.concatenate(dt2d_list)

    nsteps = vx.shape[3]

    return x, y, vx, vy, vz, elmap, dt2d, metadata, nsteps

def read_binary_file(filename, only_mesh = False, debug = False):
    # Open the file in binary mode
    print(f'Reading {filename:s}:')
    with open(filename, 'rb') as f:
        # Step 1: Read the header (assuming a fixed size of 128 bytes for example)
        header = f.read(116).split()
        if debug:
            print(header)

        # extract version and file type
        head = header[0].decode()
        if head[3] == 'f':
            if_half = False
            version = int(head[1])
        elif head[3] == 'h':
            if_half = True
            version = int(head[1])
        elif head[1:4] == 'tor':
            if_half = False
            version = 0
        else:
            print('Header filetype inconsistent.')
            print(head)
            sys.exit()
        
        # extract word size
        wdsize = int(header[1])

        # identify endian encoding
        etagb = f.read(4)
        etagL = struct.unpack('<f', etagb)[0]; etagL = int(etagL*1e5)/1e5
        etagB = struct.unpack('>f', etagb)[0]; etagB = int(etagB*1e5)/1e5
        if (etagL == 6.54321):
           emode = '<'
        elif (etagB == 6.54321):
           emode = '>'

        print(f'  read metadata')
        metadata = {
            'version': version,
            'if_half': if_half,
            'wdsize': wdsize,
            'lx1': read_int(f, emode, 1)[0],
            'ly1': read_int(f, emode, 1)[0],
            'nelf': read_int(f, emode, 1)[0],
            'time': read_flt(f, emode, wdsize, 1)[0],
            'nsave': read_int(f, emode, 1)[0],
            'lbuf': read_int(f, emode, 1)[0],
        }
        lx1 = metadata['lx1']
        ly1 = metadata['ly1']
        nsave = metadata['nsave']
        nelf = metadata['nelf']

        print(f'  read elmap')
        elmap = read_int(f,emode,nelf)
        print(f'  read dt')
        dt2d = read_flt(f,emode,wdsize,nsave)
        idx = np.argsort(elmap)
        
        nxy = lx1 * ly1
        if debug:
            xavg, yavg = np.zeros(nelf,), np.zeros(nelf,)
        print(f'  read x')
        x = np.zeros((lx1,ly1,nelf))
        for i in idx:
            xel = read_flt(f,emode,wdsize,nxy)
            x[:,:,i] = xel.reshape((lx1,ly1), order='F')
            if debug:
                xavg[i] = np.mean(xel)
        print(f'  read y')
        y = np.zeros((lx1,ly1,nelf))
        for i in idx:
            yel = read_flt(f,emode,wdsize,nxy)
            y[:,:,i] = yel.reshape((lx1,ly1), order='F')
            if debug:
                yavg[i] = np.mean(yel)
        vx, vy, vz = np.empty((lx1,ly1,nelf,nsave)), np.empty((lx1,ly1,nelf,nsave)), np.empty((lx1,ly1,nelf,nsave))
        if not only_mesh:
            print(f'  read vxyz for {nsave:d} snapshots.')
            if debug:
                print(f'  Read fld:')
            for ibuf in range(nsave):
                if debug:
                    if ibuf % 20 == 0:
                        print(f'\n  ', end='')
                    print(f' {ibuf+1:3d}', end='')
                for i in idx:
                    eldata = read_flt(f,emode,wdsize,nxy)
                    vx[:,:,i,ibuf] = eldata.reshape((lx1,ly1), order='F')
                for i in idx:
                    eldata = read_flt(f,emode,wdsize,nxy)
                    vy[:,:,i,ibuf] = eldata.reshape((lx1,ly1), order='F')
                for i in idx:
                    eldata = read_flt(f,emode,wdsize,nxy)
                    vz[:,:,i,ibuf] = eldata.reshape((lx1,ly1), order='F')
            if debug:
                print(f'')

        # Return the parsed data
        return x, y, vx, vy, vz, elmap, dt2d, metadata