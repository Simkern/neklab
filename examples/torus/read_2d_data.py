import sys, os
import struct, glob
import numpy as np

class data2d:

    def __init__(self, x, y, vx, vy, vz, dt, elmap, if_sol=True):
        self.x     = x
        self.y     = y
        self.vx    = vx
        self.vy    = vy
        self.vz    = vz
        self.dt    = dt
        self.elmap = elmap
        self.if_sol= if_sol
    
    def check_dims(self, dims):
        passed = True
        if not (self.vx.shape == self.vy.shape == self.vz.shape == dims):
            print(f'Velocity arrays have inconsistent sizes.')
            print(f'vx, vy, vz:')
            print(self.vx.shape)
            print(self.vy.shape)
            print(self.vz.shape)
            print(dims)
            passed = False
        if not (len(self.dt) == dims[-1]):
            print(f'timestep array has an inconsistent size.')
            print(f'dt:')
            print(len(self.dt))
            print(dims[3])
            passed = False
        return passed

class meta2d:

    def __init__(self, version, if_half, wdsize, emode, lx1, ly1, nelf, time, nsave, lbuf):
        self.version = version
        self.if_half = if_half
        self.wdsize  = wdsize
        self.emode   = emode
        self.lx1     = lx1
        self.ly1     = ly1
        self.nelf    = nelf
        self.time    = time
        self.nsave   = nsave
        self.lbuf    = lbuf

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

    pat = os.path.join(cwd, filepattern+'[0-9][0-9][0-9].fld')

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

    vx_list, vy_list, vz_list, dt_list = [], [], [], []
    
    if only_mesh:
        data_part, meta = read_binary_file(file, only_mesh=True)
        vx = data_part.vx
        vy = data_part.vy
        vz = data_part.vz
        dt = data_part.dt
        nsteps = 0
        is_sol = False
    else:
        for file in files:
            data_part, meta = read_binary_file(file, only_mesh=False)
            # Append the velocity arrays along with other necessary checks
            vx_list.append(data_part.vx)
            vy_list.append(data_part.vy)
            vz_list.append(data_part.vz)
            dt_list.append(data_part.dt)

        # Concatenate the velocity arrays along the last dimension (axis=3)
        vx = np.concatenate(vx_list, axis=3)
        vy = np.concatenate(vy_list, axis=3)
        vz = np.concatenate(vz_list, axis=3)
        dt = np.concatenate(dt_list)
        nsteps = vx.shape[3]
        is_sol = True

    data = data2d(data_part.x, data_part.y, vx, vy, vz, dt, data_part.elmap, is_sol)

    return data, meta, nsteps

def read_binary_file(filename, only_mesh = False, debug = False):
    # Open the file in binary mode
    print(f'\nReading {filename:s}:')
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
        print('  read metadata ', end='')
        lx1   = read_int(f, emode, 1)[0]
        ly1   = read_int(f, emode, 1)[0]
        nelf  = read_int(f, emode, 1)[0]
        time  = read_flt(f, emode, wdsize, 1)[0]
        nsave = read_int(f, emode, 1)[0]
        lbuf  = read_int(f, emode, 1)[0]
        meta = meta2d(version, if_half, wdsize, emode, lx1, ly1, nelf, time, nsave, lbuf)

        print('elmap ', end='')
        elmap = read_int(f,emode,nelf)
        print('dt')
        dt2d = read_flt(f,emode,wdsize,nsave)
        idx = np.argsort(elmap)
        
        nxy = lx1 * ly1
        print('  read x ', end='')
        x = np.zeros((lx1,ly1,nelf))
        for i in idx:
            xel = read_flt(f,emode,wdsize,nxy)
            x[:,:,i] = xel.reshape((lx1,ly1), order='F')
        print('y')
        y = np.zeros((lx1,ly1,nelf))
        for i in idx:
            yel = read_flt(f,emode,wdsize,nxy)
            y[:,:,i] = yel.reshape((lx1,ly1), order='F')
        vx, vy, vz = np.empty((lx1,ly1,nelf,nsave)), np.empty((lx1,ly1,nelf,nsave)), np.empty((lx1,ly1,nelf,nsave))
        if_sol = False
        if not only_mesh:
            print(f'  read vxyz for {nsave:d} snapshot(s).')
            if_sol = True
            if debug:
                print('  Read fld:')
            for ibuf in range(nsave):
                if debug:
                    if ibuf % 20 == 0:
                        print('\n  ', end='')
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
                print('')
        data = data2d(x, y, vx, vy, vz, dt2d, elmap, if_sol)
        print('')

        # Return the parsed data
        return data, meta