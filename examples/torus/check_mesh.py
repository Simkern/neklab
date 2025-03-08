import os, sys, glob, copy
import numpy as np
import matplotlib.pyplot as plt
import pymech as pm

from read_2d_data import read_binary_file
from plot_2d_data import plot_2d_fld
from manipulate_2d_data import symmetrize_fld, get_ord, flip_data

def measure_symmetry_error(re2name):
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
            #if mesh.ndim > 2:
            #   print(f'Element {iel+1:3d} (left) : {xcl[il]:8.4f} {ycl[il]:8.4f} {zc:8.4f}')
            #else:
            #   print(f'Element {iel+1:3d} (left) : {xcl[il]:8.4f} {ycl[il]:8.4f}')
            il += 1
         else:
            xcr[ir] = xc
            ycr[ir] = yc
            elidr.append(iel)
            #if mesh.ndim > 2:
            #   print(f'Element {iel+1:3d} (right): {xcr[ir]:8.4f} {ycr[ir]:8.4f} {zc:8.4f}')
            #else:
            #   print(f'Element {iel+1:3d} (right): {xcr[ir]:8.4f} {ycr[ir]:8.4f}')
            ir += 1

   l = xcl + ycl
   r = -xcr + ycr
   found = [ False for i in range(int(mesh.nel/2)) ]
   errmax = 0.0

   toplot = []

   for i, il in enumerate(l):
      #print(f'Element {i+1:3d}:')
      for j, jr in enumerate(r):
         yerr = abs(ycl[i] - ycr[j])
         err  = abs(il-jr)
         if yerr < 1e-6 and err < 1e-6:
            if not found[i]:
               if err > 1e-5:
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
               print(f'\tElement {j+1:3d}:')
               print(f'\t  c0: {xcl[i]:16.12f} {ycl[i]:16.12f}')
               print(f'\t  c1: {xcr[j]:16.12f} {ycr[j]:16.12f}')
               print(f'\t  error {abs(il-jr):16.12e}')
               print('Second element found!')
               sys.exit()

   print(f'Maximum error: {errmax:16.12e}')

if __name__ == "__main__":

   re2name = 'geom/test2D.re2'
   measure_symmetry_error(re2name)

   dir = 'steady/newton'
   ddirs = sorted([ d for d in glob.glob(os.path.join(dir,'v*')) if os.path.isdir(d) ])

   # symmetry error
   basename = 'torus.re2'
   for d in ddirs:
      if 'full' in d:
         re2name = os.path.join(d, basename)
         print(f'\n{d}: ', end='')
         measure_symmetry_error(re2name)

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
         mesh = d.split('_')[0]
         fld_s, fld_a = symmetrize_fld(data, data.vx)
         fig, ax = plt.subplots(2, 1)
         plot_2d_fld(ax[0], data.x, data.y, fld_s, draw_elements=True, draw_mesh=True)
         ax[0].set_title(f'Symmetric part: min/max = {fld_s.min():10.2e}/{fld_s.max():10.2e}')
         plot_2d_fld(ax[1], data.x, data.y, fld_a, draw_elements=True, draw_mesh=True)
         ax[1].set_title(f'Antisymmetric part: min/max = {fld_a.min():10.2e}/{fld_a.max():10.2e}')
         fig.suptitle(f'Mesh {mesh}: PO {dict['PO']}') 
   #print(datav)

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

   # Correct loss of precision in extruded 3D mesh for the half pipe
   nslice   = mesh2Dh.nel
   nslice3D = [ iel % nslice for iel in range(mesh3Dh.nel) ]
   for iel in range(mesh3Dh.nel):
      islice = nslice3D[ iel ]
      xy3D = mesh3Dh.elem[iel   ].pos[:2,0,:,:]
      xy2D = mesh2Dh.elem[islice].pos[:2,0,:,:]
      mesh3Dh.elem[iel].pos[:2,0,:,:] = xy2D
      mesh3Dh.elem[iel].pos[:2,1,:,:] = xy2D
      #print(f'Element {iel:3d}: err = {np.sum(abs(xy3D-xy2D)):8.2e}')

   # Enforce symmetry on the 2D full pipe mesh
   # get centroid location of elements in the half pipe mesh
   rh, phih = np.empty((mesh2Dh.nel,)), np.empty((mesh2Dh.nel,))
   for iel, elh in enumerate(mesh2Dh.elem):
      xc, yc, _ = elh.centroid
      rh[iel]   = np.sqrt(xc**2 + yc**2)
      phih[iel] = np.atan2(yc, xc)

   f2h = np.array([ -1 for i in range(mesh2Df.nel) ], dtype=int)
   # find corresponding mapping to full pipe mesh
   for ielf, elf in enumerate(mesh2Df.elem):
      xc, yc, zc = elf.centroid
      s    = xc/abs(xc) # flip the elements on the left half plane
      rf   = np.sqrt(xc**2 + yc**2)
      phif = np.atan2(yc, s*xc)
      #print(f'Element {iel:3d}:')
      for ielh, (r,phi,elh) in enumerate(zip(rh,phih,mesh2Dh.elem)):
         err = abs(rf - r) + abs(phif - phi)
         #print(f'\t{ielh:3d}: err = {err:8.2e}')
         if err < 1e-6:
            f2h[ielf] = ielh
            x0 = s*elf.pos[0,0,:,:].squeeze()
            y0 =   elf.pos[1,0,:,:].squeeze()
            x1 =   elh.pos[0,0,:,:].squeeze()
            y1 =   elh.pos[1,0,:,:].squeeze()
            # map elements, use loose tolerance!
            order = get_ord(x0,y0,x1,y1,tol=1e-4)
            # update the x,y data of the full pipe
            mesh2Df.elem[ielf].pos[0,0,:,:] = s*flip_data(x1,order)
            mesh2Df.elem[ielf].pos[1,0,:,:] =   flip_data(y1,order)
            break

   if f2h.min() < 0:
      print(f'Error: Not all elements mapped.')
      sys.exit()
   
   # update the 3D mesh based on the 2D one for the full pipe
   nslice   = mesh2Df.nel
   nslice3D = [ iel % nslice for iel in range(mesh3Df.nel) ]
   for iel in range(mesh3Df.nel):
      islice = nslice3D[ iel ]
      xy3D = mesh3Df.elem[iel   ].pos[:2,0,:,:]
      xy2D = mesh2Df.elem[islice].pos[:2,0,:,:]
      mesh3Df.elem[iel].pos[:2,0,:,:] = xy2D
      mesh3Df.elem[iel].pos[:2,1,:,:] = xy2D
      #print(f'Element {iel:3d}: err = {np.sum(abs(xy3D-xy2D)):8.2e}')

   ax = plt.figure().add_subplot(projection='3d')
   for el in mesh3Dh.elem:
      xc, yc, zc = el.centroid
      ax.scatter(xc,yc,zc,s=30,c=r)
   ax = plt.figure().add_subplot(projection='3d')
   for el in mesh3Df.elem:
      xc, yc, zc = el.centroid
      ax.scatter(xc,yc,zc,s=30,c=r)
   fig, axs = plt.subplots(1,2)
   ax = axs[0]
   for el in mesh2Dh.elem:
      xc, yc, _ = el.centroid
      ax.scatter(xc,yc,s=30,c=r)
   ax = axs[1]
   for el in mesh2Df.elem:
      xc, yc, _ = el.centroid
      ax.scatter(xc,yc,s=30,c=r)
   plt.show()
   
      
