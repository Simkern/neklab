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
    with open(file_path, 'rb') as f:
        # Step 1: Read the header (assuming a fixed size of 128 bytes for example)
        header = f.read(116).split()

        # extract word size
        wdsize = int(header[1])

        print(header)

        # identify endian encoding
        etagb = f.read(4)
        etagL = struct.unpack('<f', etagb)[0]; etagL = int(etagL*1e5)/1e5
        etagB = struct.unpack('>f', etagb)[0]; etagB = int(etagB*1e5)/1e5
        if (etagL == 6.54321):
           emode = '<'
        elif (etagB == 6.54321):
           emode = '>'

        lx1 = read_int(f,emode,1)[0]
        ly1 = read_int(f,emode,1)[0]
        nelf = read_int(f,emode,1)[0]
        time = read_flt(f,emode,wdsize,1)[0]
        nsave = read_int(f,emode,1)[0]
        lbuf = read_int(f,emode,1)[0]

        elmap = read_int(f,emode,nelf)
        # Step 5: Read the x coordinates (lx1 * ly1 * nelf reals)
        npts = lx1 * ly1 * nelf
        x  = read_flt(f,emode,wdsize,npts)
        y  = read_flt(f,emode,wdsize,npts)
        vx, vy, vz = np.empty((lx1*ly1*nelf, nsave)), np.empty((lx1*ly1*nelf, nsave)), np.empty((lx1*ly1*nelf, nsave))
        for i in range(nsave):
            vx[:,i] = read_flt(f,emode,wdsize,npts)
            vy[:,i] = read_flt(f,emode,wdsize,npts)
            vz[:,i] = read_flt(f,emode,wdsize,npts)

        # Return the parsed data
        return header, emode, lx1, ly1, nelf, x, y, vx, vy, vz

# Example usage:
file_path = '01run/c000.fld'
header, endian_test_string, lx1, ly1, nelf, x, y, vx, vy, vz = read_binary_file(file_path)

plt.figure()
plt.scatter(x,y,50,'k')

plt.figure()
plt.plot(x)
plt.plot(y)

plt.figure()
x2d = np.reshape(x,  [lx1,ly1,nelf], order='F')
y2d = np.reshape(y,  [lx1,ly1,nelf], order='F')
v2d = np.reshape(vx, [lx1,ly1,nelf], order='F')
for i in range(nelf):
    plt.pcolor(x2d[:,:,i],y2d[:,:,i],v2d[:,:,i])
plt.colorbar()
plt.show()