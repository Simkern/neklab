import struct
import numpy as np
import matplotlib.pyplot as plt

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

def read_binary_file(file_path):
    # Open the file in binary mode
    print(f'Reading {file_path:s}:')
    with open(file_path, 'rb') as f:
        # Step 1: Read the header (assuming a fixed size of 128 bytes for example)
        header = f.read(116).split()

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

        print(f'\tread metadata')
        lx1 = read_int(f,emode,1)[0]
        ly1 = read_int(f,emode,1)[0]
        nelf = read_int(f,emode,1)[0]
        time = read_flt(f,emode,wdsize,1)[0]
        nsave = read_int(f,emode,1)[0]
        lbuf = read_int(f,emode,1)[0]

        print(f'\tread elmap')
        elmap = read_int(f,emode,nelf)
        print(f'\tdt information')
        dt2d = read_flt(f,emode,wdsize,nsave)
        idx = np.argsort(elmap)
        
        #print(" ".join([f'{elmap[i]:d}' for i in idx]))
        nxy = lx1*ly1
        print(f'\tread x')
        x = np.zeros((lx1,ly1,nelf))
        xavg = np.zeros(nelf,)
        for i in idx:
            xel = read_flt(f,emode,wdsize,nxy)
            x[:,:,i] = xel.reshape((lx1,ly1), order='F')
            xavg[i] = np.mean(xel)
        print(f'\tread y')
        y = np.zeros((lx1,ly1,nelf))
        yavg = np.zeros(nelf,)
        for i in idx:
            yel = read_flt(f,emode,wdsize,nxy)
            y[:,:,i] = yel.reshape((lx1,ly1), order='F')
            yavg[i] = np.mean(yel)
        vx, vy, vz = np.empty((lx1,ly1,nelf,nsave)), np.empty((lx1,ly1,nelf,nsave)), np.empty((lx1,ly1,nelf,nsave))
        print(f'\tread vxyz for {nsave:d} snapshots.')
        for ibuf in range(nsave):
            #print(f'\tread vx')
            for id, i in enumerate(idx):
                eldata = read_flt(f,emode,wdsize,nxy)
                vx[:,:,i,ibuf] = eldata.reshape((lx1,ly1), order='F')
            #print(f'\tread vy')
            for id, i in enumerate(idx):
                eldata = read_flt(f,emode,wdsize,nxy)
                vy[:,:,i,ibuf] = eldata.reshape((lx1,ly1), order='F')
            #print(f'\tread vz')
            for id, i in enumerate(idx):
                eldata = read_flt(f,emode,wdsize,nxy)
                vz[:,:,i,ibuf] = eldata.reshape((lx1,ly1), order='F')

        # Return the parsed data
        return header, emode, lx1, ly1, nelf, elmap, x, y, vx, vy, vz