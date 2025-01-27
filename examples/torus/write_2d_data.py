import sys, os
import struct
import glob
import numpy as np

def write_int(f, emode, nvar, llist):
    """Write integer array"""
    packed_data = struct.pack(emode + nvar * 'i', *llist)
    f.write(packed_data)

def write_flt(f, emode, wdsize, nvar, llist):
    """Write real (floating-point) array"""
    if wdsize == 4:
        realtype = 'f'  # 4-byte float
    elif wdsize == 8:
        realtype = 'd'  # 8-byte double
    packed_data = np.array(llist, dtype=emode + realtype)
    packed_data.tofile(f)  # Write the array directly to the file

def write_fields(filepattern, x, y, vx, vy, vz, elmap, dt2d, metadata, nsteps, cwd='.'):

    lx1 = metadata['lx1']
    ly1 = metadata['ly1']
    nelf = metadata['nelf']
    lbuf = metadata['lbuf']
    # consistency check
    nsteps = len(dt2d)
    dims = [ lx1, ly1, nelf, nsteps ]
    if not (vx.shape == vy.shape == vz.shape == dims):
        print(f'Velocity arrays have inconsistent sizes.')
        print(f'vx, vy, vz, dt2d:')
        print(vx.shape)
        print(vy.shape)
        print(vz.shape)
        print(dt2d.shape)
        sys.exit()

    if nsteps > lbuf:
        vxpart = [ vx[..., i:i+lbuf] for i in range(0, nsteps, lbuf) ]
        vypart = [ vy[..., i:i+lbuf] for i in range(0, nsteps, lbuf) ]
        vzpart = [ vz[..., i:i+lbuf] for i in range(0, nsteps, lbuf) ]
        dtpart = [ dt2d[i:i+lbuf] for i in range(0, nsteps, lbuf) ]
        nsave_list = [ len(d) for d in dtpart ]
    else:
        vxpart = vx
        vzpart = vy
        vzpart = vz
        dtpart = dt2d
        nsave_list = nsteps

    if not os.path.isdir(cwd):
        print(f"The path '{cwd}' does not exist. Abort.")
        sys.exit()
    
    files = glob.glob(os.path.join(cwd,filepattern+'*.fld'))

    if len(files) > 0:
        print('The following files exist and will be overwritten:')
        for file in files:
            print(f'   {file}')
        user_input = input(f"OK? (yes/no) [y] : ")
        if user_input.lower() not in ['', 'y', 'yes']:
            sys.exit()

    for ifile, (vxp, vyp, vzp, dtp, nsave) in enumerate(zip(vxpart, vypart, vzpart, dtpart, nsave_list)):
        filename = os.path.join(cwd, filepattern+f'{ifile+1:3d}'+'.fld')
        metadata['nsave'] = nsave
        write_binary_file(filename, x, y, vxp, vyp, vzp, elmap, dtp, metadata)

def write_binary_file(filename, x, y, vx, vy, vz, elmap, dt2d, metadata, debug=False):

    version = metadata['version']
    if_half = metadata['if_half']
    wdsize = metadata['wdsize']
    lx1 = metadata['lx1']
    ly1 = metadata['ly1']
    nelf = metadata['nelf']
    time = metadata['time']
    nsave = metadata['nsave']
    lbuf = metadata['lbuf']
    emode = metadata['emode']
    nxy = lx1*ly1

    dims = [ lx1, ly1, nelf, nsave ]
    if not (vx.shape == vy.shape == vz.shape == dims) or not len(dt2d) == nsave:
        print(f'Velocity arrays have inconsistent sizes.')
        print(f'vx, vy, vz, dt2d:')
        print(vx.shape)
        print(vy.shape)
        print(vz.shape)
        print(dt2d.shape)
        sys.exit()
    
    with open(filename, 'wb') as f:
        # Step 1: Write the header (116 bytes)
        if if_half:
            head = f'#{version}th'.encode()
        else:
            head = f'#{version}tf'.encode()
        
        # Write the header data
        f.write(head)
        
        # Write word size and endianness data (emulating reading endian encoding from the header)
        write_int(f, emode, 3, [lx1, ly1, nelf])
        write_flt(f, emode, wdsize, 1, time)
        write_int(f, nsave, 1, lbuf)

        # Write element mapping and dt2d data
        write_int(f, emode, nelf, elmap)
        write_flt(f, emode, wdsize, nsave, dt2d)

        if debug:
            print(f"Written metadata: {metadata}")
            print(f"Written elmap: {elmap}")
            print(f"Written dt2d: {dt2d}")
        
        print(f'  write x')
        for i in elmap:
            write_flt(f, emode, wdsize, nxy, x[:,:,i].flatten(order='F'))
            
        print(f'  write y')
        for i in elmap:
            write_flt(f, emode, wdsize, nxy, y[:,:,i].flatten(order='F'))

        # Write vx, vy, vz for the snapshots
        print(f'  write vxyz for {nsave:d} snapshots.')
        for ibuf in range(nsave):
            for i in elmap:
                write_flt(f, emode, wdsize, nxy, vx[:,:,i,ibuf].flatten(order='F'))
            for i in elmap:
                write_flt(f, emode, wdsize, nxy, vy[:,:,i,ibuf].flatten(order='F'))
            for i in elmap:
                write_flt(f, emode, wdsize, nxy, vz[:,:,i,ibuf].flatten(order='F'))