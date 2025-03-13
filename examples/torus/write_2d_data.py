import sys, os
import struct
import glob
import numpy as np
from read_2d_data import data2d, meta2d

def write_int(f, emode, nvar, llist):
    """Write integer array"""
    packed_data = struct.pack(emode + nvar*'i', *llist)
    f.write(packed_data)

def write_flt(f, emode, wdsize, nvar, llist):
    """Write real (floating-point) array"""
    if wdsize == 4:
        realtype = 'f'  # 4-byte float
    elif wdsize == 8:
        realtype = 'd'  # 8-byte double
    packed_data = struct.pack(emode+nvar*realtype, *llist)
    f.write(packed_data)

def write_fields(filepattern, data, meta, nsteps, cwd='.', force=False):

    lx1  = meta.lx1
    ly1  = meta.ly1
    nelf = meta.nelf
    lbuf = meta.lbuf
    # consistency check
    nsteps = len(data.dt)
    dims = ( lx1, ly1, nelf, nsteps )
    
    if not data.check_dims(dims):
        print('Error in write_fields.')
        sys.exit()

    if nsteps > lbuf:
        vxpart = [ data.vx[..., i:i+lbuf] for i in range(0, nsteps, lbuf) ]
        vypart = [ data.vy[..., i:i+lbuf] for i in range(0, nsteps, lbuf) ]
        vzpart = [ data.vz[..., i:i+lbuf] for i in range(0, nsteps, lbuf) ]
        dtpart = [ data.dt[i:i+lbuf] for i in range(0, nsteps, lbuf) ]
        nsave_list = [ len(d) for d in dtpart ]
        nfiles = len(dtpart)
    else:
        vxpart = [ data.vx ]
        vypart = [ data.vy ]
        vzpart = [ data.vz ]
        dtpart = [ data.dt ]
        nsave_list = [ nsteps ]
        nfiles = 1

    if not os.path.isdir(cwd):
        print(f"The path '{cwd}' does not exist. Abort.")
        sys.exit()
    
    files = glob.glob(os.path.join(cwd,filepattern+'*.fld'))

    if len(files) > 0 and not force:
        print('The following files exist and will be overwritten:')
        for file in files:
            print(f'   {file}')
        user_input = input(f"OK? (yes/no) [y] : ")
        if user_input.lower() not in ['', 'y', 'yes']:
            sys.exit()

    for ifile in range(nfiles):
        filename = os.path.join(cwd, filepattern+f'{ifile+1:03d}'+'.fld')
        dsave = data2d(data.x, data.y, 
                   vxpart[ifile], vypart[ifile], vzpart[ifile], dtpart[ifile], 
                   data.elmap)
        meta.nsave = nsave_list[ifile]
        write_binary_file(filename, dsave, meta)
        print(f'{filename} written.')

def write_binary_file(filename, data, meta, debug=False):

    version = meta.version
    if_half = meta.if_half
    wdsize  = meta.wdsize
    lx1     = meta.lx1
    ly1     = meta.ly1
    nelf    = meta.nelf
    time    = meta.time
    nsave   = meta.nsave
    lbuf    = meta.lbuf
    emode   = meta.emode
    nxy     = lx1*ly1

    if not data.vx.shape == data.vx.shape == data.vx.shape:
        print('Velocity arrays have inconsistent sizes.')
        print('vx, vy, vz:')
        print(data.vx.shape)
        print(data.vy.shape)
        print(data.vz.shape)
        sys.exit()

    if nsave > lbuf:
        print(f'nsave = {nsave} > {lbuf} = lbuf. Abort.')
        sys.exit()
    elif nsave > data.vx.shape[-1]:
        print(f'nsave = {nsave} > {data.vx.shape[-1]} = vx.shape[-1].')
        nsave = data.vx.shape[-1]
        meta.nsave = nsave
        print(f'Reset nsave = {nsave}')
    elif nsave < data.vx.shape[-1]:
        print(f'nsave = {nsave} < {data.vx.shape[-1]} = vx.shape[-1].')
        print('Not all data will be written to file.')
    
    print(f'\nWriting {filename:s}:')
    with open(filename, 'wb') as f:
        # Step 1: Write the header (116 bytes)
        if if_half:
            id = 'th'
        else:
            id = 'tf'

        header = (
            f'#{version}{id} {wdsize} '
            f'(lx1, ly1 ={lx1:9d}{ly1:9d}) '
            f'(nelf ={nelf:9d}) '
            f'(time ={time:17.9e}) '
            f'(nsave, lbuf ={nsave:9d}{lbuf:9d})'
        )
        
        # Write the header datal
        f.write(header.ljust(116).encode('utf-8'))
        f.write(struct.pack(emode+'f', 6.54321))
        
        print('  write metadata ', end='')
        # Write word size and endianness data (emulating reading endian encoding from the header)
        write_int(f, emode, 3, [lx1, ly1, nelf])
        write_flt(f, emode, wdsize, 1, [time])
        write_int(f, emode, 2, [nsave, lbuf])

        # Write element mapping and dt2d data
        print('elmap ', end='')
        write_int(f, emode, nelf, data.elmap)
        print('dt')
        write_flt(f, emode, wdsize, nsave, data.dt[:nsave])

        if debug:
            print('Written metadata: meta')
            print(f'Written elmap: {data.elmap}')
            print(f'Written dt2d: {data.dt}')
        
        idx = np.argsort(data.elmap)
        print('  write x ', end='')
        for i in idx:
            write_flt(f, emode, wdsize, nxy, data.x[:,:,i].flatten(order='F'))
            
        print('y')
        for i in idx:
            write_flt(f, emode, wdsize, nxy, data.y[:,:,i].flatten(order='F'))

        # Write vx, vy, vz for the snapshots
        print(f'  write vxyz for {nsave:d} snapshot(s).')
        for ibuf in range(nsave):
            for i in idx:
                write_flt(f, emode, wdsize, nxy, data.vx[:,:,i,ibuf].flatten(order='F'))
            for i in idx:
                write_flt(f, emode, wdsize, nxy, data.vy[:,:,i,ibuf].flatten(order='F'))
            for i in idx:
                write_flt(f, emode, wdsize, nxy, data.vz[:,:,i,ibuf].flatten(order='F'))