import sys, os, re
import copy
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import pymech as pm

def add_h(s):
    return re.sub(r'(\d)(D\.rea)', r'_h\1\2', s)

def split_mesh(fname_in):
   #meshfldr = 'geom'
   meshfldr = '.'
   #fname_in  = 'torus_coarse2D.rea'
   #fname_out = 'torus_coarse2Dh.rea'
   fname_out = add_h(fname_in)
   fname = os.path.join(meshfldr,fname_in)

   print('Read mesh '+fname)
   mesh = pm.readrea(fname)

   # copy the original mesh
   hmesh = copy.deepcopy(mesh)

   # cut mesh, only keep y > 0
   # also, get centroids for plot
   helem = []
   keep = [ False for i in range(mesh.nel) ]
   xyc  = np.empty((mesh.nel,2))
   xych = np.empty((int(mesh.nel/2),2))
   for (iel, el) in enumerate(mesh.elem):
      x, y, z = el.centroid
      xyc[iel,:] = x, y
      if x > 0.0:
         keep[iel] = True
         helem.append(mesh.elem[iel])

   # update cut mesh
   hmesh.elem = helem
   hmesh.nel  = sum(keep)
   hmesh.update_ncurv()

   xych = np.empty((hmesh.nel,2))
   xcm = np.empty(0)
   ycm = np.empty(0)
   xcs = np.empty(0)
   ycs = np.empty(0)
   ms  = [ False for i in range(hmesh.nel) ]
   sym = [ False for i in range(hmesh.nel) ]
   for (iel, el) in enumerate(hmesh.elem):
      x, y, z = el.centroid
      xych[iel,:] = x, y
      has_bc = False
      for iface in range(2*hmesh.ndim):
         bc = el.bcs[0, iface][0]
         fc = el.face_center(iface)
         if bc == 'MS':
            ms[iel] = True
            #xycm[iel] = fc[:2]
            xcm = np.append(xcm, fc[0])
            ycm = np.append(ycm, fc[1])
            el.bcs[0, iface][0] = 'W  '
            has_bc = True
         if abs(fc[0]) < 1e-6:
            sym[iel] = True
            #xycs[iel] = fc[:2]
            xcs = np.append(xcs, fc[0])
            ycs = np.append(ycs, fc[1])
            el.bcs[0, iface][0] = 'SYM'
            has_bc = True
         if not has_bc:
            el.bcs[0, iface][0] = 'E'

   for iel, el in enumerate(hmesh.elem):
      print(f'El {iel+1}: ', ', '.join([ f'{c}' for c in el.centroid[:2] ]))
      for iface in range(2*hmesh.ndim):
         bc = el.bcs[0, iface][0]
         print(f'\t face {iface+1}: {bc}')

   fname = os.path.join(meshfldr,fname_out)

   print('Write mesh '+fname)
   pm.writerea(fname, hmesh)

   fig, ax = plt.subplots(1,1, figsize=(10, 8))
   ax.scatter(xyc[:,0], xyc[:,1], c = 'b', marker='o')
   ax.scatter(xych[:,0], xych[:,1], c = 'r', marker='o')
   ax.scatter(xcm, ycm, c = 'k', marker='x')
   ax.scatter(xcs, ycs, c = 'k', marker='s')

   plt.show()

if __name__ == "__main__":
   split_mesh('torus_coarse2D.rea')