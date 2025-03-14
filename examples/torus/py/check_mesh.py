import os, sys, glob, copy
import numpy as np
import matplotlib.pyplot as plt
import pymech as pm

from read_2d_data import read_binary_file
from write_2d_data import write_binary_file
from plot_2d_data import plot_2d_fld
from manipulate_2d_data import symmetrize_fld, get_ord, flip_data

def repair_extruded_mesh(m2D, m3D, plot=False, verb=False):
   m3D_new = copy.deepcopy(m3D)
   # Correct loss of precision in extruded 3D mesh for the half pipe
   nslice = m2D.nel
   nslice3D = [ iel % nslice for iel in range(m3D.nel) ]
   for iel in range(m3D.nel):
      islice = nslice3D[ iel ]
      xy2D = m2D.elem[islice].pos[:2,0,:,:]
      m3D_new.elem[iel].pos[:2,0,:,:] = xy2D
      m3D_new.elem[iel].pos[:2,1,:,:] = xy2D
      if verb:
         xy3D = m3D.elem[iel].pos[:2,0,:,:]
         print(f'Element {iel:3d}: max err = {abs(xy3D-xy2D).max():8.2e}')
      
   # Check resulting mesh
   max_coord_shift = 0.0
   if plot:
      fig = plt.figure()
      ax = fig.add_subplot(projection='3d')
   for iel, (el, elc) in enumerate(zip(m3D.elem,m3D_new.elem)):
      x , y , z  = el.centroid
      xc, yc, zc = elc.centroid
      err = abs(el.pos - elc.pos).max()
      max_coord_shift = max(err, max_coord_shift)
      if plot:
         ax.scatter(x , y , z , s=10, c='r', label='original')
         ax.scatter(xc, yc, zc, s=50, c='k', marker='+', label='updated')
      
   print(f'Maximum coordinate shift: {max_coord_shift:16.12e}')

   if plot:
      handles, labels = ax.get_legend_handles_labels()
      by_label = dict(zip(labels, handles))
      plt.legend(by_label.values(), by_label.keys())
      plt.show()
      
   return m3D_new

def enforce_symmetry_2D(m_half, m_full, plot=False, verb=False):
   # Enforce symmetry on the 2D full pipe mesh
   if m_half.ndim == 3 or m_full.ndim == 3:
      print(f'The input meshes must be 2D.')
      sys.exit()
   # get centroid location of elements in the half pipe mesh
   rh, phih = np.empty((m_half.nel,)), np.empty((m_half.nel,))
   for iel, elh in enumerate(m_half.elem):
      xc, yc, _ = elh.centroid
      rh[iel]   = np.sqrt(xc**2 + yc**2)
      phih[iel] = np.atan2(yc, xc)

   m_full_new = copy.deepcopy(m_full)
   f2h = np.array([ -1 for i in range(m_full.nel) ], dtype=int)
   # find corresponding mapping to full pipe mesh
   for ielf, elf in enumerate(m_full.elem):
      xc, yc, _ = elf.centroid
      s    = xc/abs(xc) # flip the elements on the left half plane
      rf   = np.sqrt(xc**2 + yc**2)
      phif = np.atan2(yc, s*xc)
      if verb:
         print(f'Element {iel:3d}:')
      for ielh, (r, phi, elh) in enumerate(zip(rh, phih, m_half.elem)):
         err = abs(rf - r) + abs(phif - phi)
         if err < 1e-6:
            if verb:
               print(f'\tmapped to element {ielh:3d}: err = {err:8.2e}')
            f2h[ielf] = ielh
            x0 = s*elf.pos[0,0,:,:].squeeze(); y0 = elf.pos[1,0,:,:].squeeze()
            x1 =   elh.pos[0,0,:,:].squeeze(); y1 = elh.pos[1,0,:,:].squeeze()
            # map elements, use loose tolerance!
            order = get_ord(x0,y0,x1,y1,tol=1e-4)
            # update the x,y data of the full pipe
            m_full_new.elem[ielf].pos[0,0,:,:] = s*flip_data(x1,order)
            m_full_new.elem[ielf].pos[1,0,:,:] =   flip_data(y1,order)
            break

   # Check that all elements have been mapped
   if f2h.min() < 0:
      print(f'Error: Not all elements mapped.')
      sys.exit()

   max_coord_shift = 0.0
   if plot:
      ax = plt.figure().add_subplot()
   for iel, (el, elc) in enumerate(zip(m_full.elem,m_full_new.elem)):
      x , y , _ = el.centroid
      xc, yc, _ = elc.centroid
      err = abs(el.pos - elc.pos).max()
      max_coord_shift = max(err, max_coord_shift)
      if plot:
         ax.scatter(x , y , s=10, c='r', label='original')
         ax.scatter(xc, yc, s=50, c='k', marker='+', label='updated')
      if verb:
         print(f'Element {iel:3d}: {err:16.12e}')
   
   print(f'Maximum coordinate shift: {max_coord_shift:16.12e}')

   if plot:
      handles, labels = ax.get_legend_handles_labels()
      by_label = dict(zip(labels, handles))
      plt.legend(by_label.values(), by_label.keys())
      plt.show()

   return m_full_new

def measure_symmetry_error(re2name, verb=False):
   mesh = pm.neksuite.readre2(re2name)

   # get element centers
   xcl = np.empty((int(mesh.nel/2),))
   xcr = np.empty((int(mesh.nel/2),))
   ycl = np.empty((int(mesh.nel/2),))
   ycr = np.empty((int(mesh.nel/2),))
   elidl = []
   elidr = []
   il = 0
   ir = 0
   for iel, el in enumerate(mesh.elem):
      xc, yc, zc = el.centroid
      flag = True
      if mesh.ndim > 2:
         zmin = 1.0
         zmax = 0.0
         for i in range(2*mesh.ndim):
            _, _, zfc = el.face_center(i)
            zmin = min(zmin, zfc)
            zmax = max(zmax, zfc)
         if zmin > 0.0:
            flag = False
      if flag:
         if xc < 0.0:
            xcl[il] = xc
            ycl[il] = yc
            elidl.append(iel)
            if verb:
               if mesh.ndim > 2:
                  print(f'Element {iel+1:3d} (left) : {xcl[il]:8.4f} {ycl[il]:8.4f} {zc:8.4f}')
               else:
                  print(f'Element {iel+1:3d} (left) : {xcl[il]:8.4f} {ycl[il]:8.4f}')
            il += 1
         else:
            xcr[ir] = xc
            ycr[ir] = yc
            elidr.append(iel)
            if verb:
               if mesh.ndim > 2:
                  print(f'Element {iel+1:3d} (right): {xcr[ir]:8.4f} {ycr[ir]:8.4f} {zc:8.4f}')
               else:
                  print(f'Element {iel+1:3d} (right): {xcr[ir]:8.4f} {ycr[ir]:8.4f}')
            ir += 1

   l = xcl + ycl
   r = -xcr + ycr
   found = [ False for i in range(int(mesh.nel/2)) ]
   errmax = 0.0

   toplot = []

   for i, il in enumerate(l):
      if verb:
         print(f'Element {i+1:3d}:')
      for j, jr in enumerate(r):
         yerr = abs(ycl[i] - ycr[j])
         err  = abs(il-jr)
         if yerr < 1e-6 and err < 1e-6:
            if not found[i]:
               if verb:
                  toplot.append(elidl[i])
                  toplot.append(elidr[j])
                  print(f'\tElement {j+1:3d}:')
                  print(f'\t  c0: {xcl[i]:16.12f} {ycl[i]:16.12f}')
                  print(f'\t  c1: {xcr[j]:16.12f} {ycr[j]:16.12f}')
                  print(f'\t  r:  {np.sqrt(xcr[j]**2+ycr[j]**2):16.12f}')
                  print(f'\t  error {abs(il-jr):16.12e}')
               found[i] = True
               errmax = max(errmax, err)
            else:
               print('Error: Second element match found!')
               print(f'\tElement {j+1:3d}:')
               print(f'\t  c0: {xcl[i]:16.12f} {ycl[i]:16.12f}')
               print(f'\t  c1: {xcr[j]:16.12f} {ycr[j]:16.12f}')
               print(f'\t  error {abs(il-jr):16.12e}')
               sys.exit()

   print(f'Maximum error: {errmax:16.12e}')
   return errmax

if __name__ == "__main__":

   #re2name = 'geom/test2D.re2'
   #errmax = measure_symmetry_error(re2name)

   dir = 'steady/newton'
   ddirs = sorted([ d for d in glob.glob(os.path.join(dir,'v*')) if os.path.isdir(d) and os.path.isfile(os.path.join(d, 'c2dtorus001.fld'))])
   outdir = 'c2d_steady'

   # symmetry error
   basename = 'torus.re2'
   for d in ddirs:
      if 'full' in d:
         re2name = os.path.join(d, basename)
         print(f'\n{d}: ', end='')
         errmax = measure_symmetry_error(re2name)

   fig, axs = plt.subplots(3,2)
   datav = {}
   for d in ddirs:
      print(d)
      fname = os.path.join(d, 'c2dtorus001.fld')
      
      data, meta = read_binary_file(fname)

      icol = 0
      full = False
      if 'full' in d:
         icol = 1
         full = True
      if 'P08' in d:
         irow = 0
         PO = 8
      elif 'P10' in d:
         irow = 1
         PO = 10
      elif 'P12' in d:
         irow = 2
         PO = 12
      else:
         print('Error.')
         sys.exit()
      ax = axs[irow,icol]
      dict = {
         'data': data,
         'meta': meta,
         'PO': PO,
         'full': full
      }
      datav[d] = dict
      if 'v0' in d:
         plot_2d_fld(ax, data.x, data.y, data.vx, draw_elements=True, draw_mesh=True)

   for d in ddirs:
      dict = datav[d]
      if dict['full']:
         data = dict['data']
         meta = dict['meta']
         mesh = d.split('_')[0].strip(dir)
         vx_s, vx_a = symmetrize_fld(data, data.vx)
         vy_s, vy_a = symmetrize_fld(data, data.vy)
         vz_s, vz_a = symmetrize_fld(data, data.vz)
         data.vx = vx_s
         data.vy = vy_s
         data.vz = vz_s
         fig, ax = plt.subplots(2, 1)
         plot_2d_fld(ax[0], data.x, data.y, vx_s, draw_elements=True, draw_mesh=False)
         ax[0].set_title(f'Symmetric part: min/max = {vx_s.min():10.2e}/{vx_s.max():10.2e}')
         plot_2d_fld(ax[1], data.x, data.y, vx_a, draw_elements=True, draw_mesh=False)
         ax[1].set_title(f'Antisymmetric part: min/max = {vx_a.min():10.2e}/{vx_a.max():10.2e}')
         fig.suptitle(f'Mesh {mesh}: PO {dict['PO']}')
         # save symm files
         fname = os.path.join(outdir, 'c2dtorus_'+mesh+f'_f_P{dict['PO']:02d}_sym.fld')
         if not os.path.isfile(fname):
            write_binary_file(fname, data, meta)
   #print(datav)
   plt.show()
   sys.exit()

   mesh2Dh_re2file = os.path.join('geomh','torus_v2_2D.re2')
   print(f'Read {mesh2Dh_re2file}')
   mesh2Dh = pm.neksuite.readre2(mesh2Dh_re2file)
   mesh3Dh_re2file = os.path.join('geomh','torus_v2_3z_3D.re2')
   print(f'Read {mesh3Dh_re2file}')
   mesh3Dh = pm.neksuite.readre2(mesh3Dh_re2file)
   mesh2Df_re2file = os.path.join('geom','torus_v2_2D.re2')
   print(f'Read {mesh2Df_re2file}')
   mesh2Df = pm.neksuite.readre2(mesh2Df_re2file)
   mesh3Df_re2file = os.path.join('geom','torus_v2_3z_3D.re2')
   print(f'Read {mesh3Df_re2file}')
   mesh3Df = pm.neksuite.readre2(mesh3Df_re2file)

   print('\nRepair extruded mesh for the half pipe:')
   m3Dh = repair_extruded_mesh(mesh2Dh, mesh3Dh, plot=True)
   
   print('\nEnfore symmetry on the 2D full pipe mesh:')
   m2Df = enforce_symmetry_2D(mesh2Dh, mesh2Df, plot=True)

   print('\nRepair extruded mesh for the full pipe after symmetry is enforced:')
   m3Df = repair_extruded_mesh(mesh2Df, mesh3Df, plot=True)
   
      
