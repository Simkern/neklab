import struct
import sys
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
        dt2d = read_flt(f,emode,wdsize,nsave)

        idx = np.argsort(elmap)
        
        print(" ".join([f'{elmap[i]:d}' for i in idx]))
        nxy = lx1*ly1
        x = np.zeros((lx1,ly1,nelf))
        for i in idx:
            xel = read_flt(f,emode,wdsize,nxy)
            x[:,:,i] = xel.reshape((lx1,ly1), order='F')
        y = np.zeros((lx1,ly1,nelf))
        for i in idx:
            yel = read_flt(f,emode,wdsize,nxy)
            y[:,:,i] = yel.reshape((lx1,ly1), order='F')
        vx, vy, vz = np.empty((lx1,ly1,nelf,nsave)), np.empty((lx1,ly1,nelf,nsave)), np.empty((lx1,ly1,nelf,nsave))
        for ibuf in range(nsave):
            for i in idx:
                eldata = read_flt(f,emode,wdsize,nxy)
                vx[:,:,i,ibuf] = eldata.reshape((lx1,ly1), order='F')
            for i in idx:
                eldata = read_flt(f,emode,wdsize,nxy)
                vy[:,:,i,ibuf] = eldata.reshape((lx1,ly1), order='F')
            for i in idx:
                eldata = read_flt(f,emode,wdsize,nxy)
                vz[:,:,i,ibuf] = eldata.reshape((lx1,ly1), order='F')

        # Return the parsed data
        return header, emode, lx1, ly1, nelf, elmap, x, y, vx, vy, vz

# Example usage:
file_path = '01run/c2dtorus000.fld'
header, endian_test_string, lx1, ly1, nelf, elmap, x, y, vx, vy, vz = read_binary_file(file_path)

ir = 4
ic = int(nelf/ir)
for i in range(ir):
    print(" ".join([ f'{el:d}' for el in elmap[i*ic:(i+1)*ic]]))

nxy = lx1*ly1
idx = 1

print('x2d')
print(f'   el {idx:3d}:')
for i in range(lx1):
    print(" ".join([f'{i:3d}:'] + [f'{dat:9.6f}' for dat in x[:,i,idx]]))
print('y2d')
print(f'   el {idx:3d}:')
for i in range(lx1):
    print(" ".join([f'{i:3d}:'] + [f'{dat:9.6f}' for dat in y[:,i,idx]]))
print('vx2d')
print(f'   el {idx:3d}:')
for i in range(lx1):
    print(" ".join([f'{i:3d}:'] + [f'{dat:9.6f}' for dat in vx[:,i,idx,0]]))

plt.figure()
plt.scatter(x.ravel(),y.ravel(),10,'k')
cmap = plt.cm.viridis
colors = [cmap(i / nelf) for i in range(nelf)]
for i in range(nelf):
    plt.scatter(x[:,:,i].ravel(),y[:,:,i].ravel(),20,color=colors[i])
    #plt.scatter(x[:,:,i].ravel(),y[:,:,i].ravel(),20,c=vx[:,:,i].ravel())
    plt.plot(x[:,0,i],y[:,0,i],color='r')

plt.colorbar()

plt.figure()
plt.plot(np.mean(x,axis = (0,1)).ravel())
plt.plot(np.mean(y,axis = (0,1)).ravel())

plt.figure()
plt.plot(np.mean(vx,axis = (0,1)).ravel())
plt.plot(np.mean(vy,axis = (0,1)).ravel())
plt.plot(np.mean(vz,axis = (0,1)).ravel())
#plt.scatter(x,y,50,'k')

#plt.figure()
#plt.plot(x)
#plt.plot(y)

centers = np.ones((lx1,ly1), dtype=bool)
centers[1:-1,1:-1] = False
idx = centers.ravel()

plt.figure()
ax = plt.gca()
vmin=vx.min()
vmax=vx.max()
for i in range(nelf):
    xi = np.squeeze(x[:,:,i])
    yi = np.squeeze(y[:,:,i])
    vxi = np.squeeze(vx[:,:,i])
    plt.contourf(xi,yi,vxi, vmin=vmin, vmax=vmax)
plt.colorbar()
plt.show()