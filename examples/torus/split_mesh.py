import sys, os, re
import copy
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import pymech as pm
import pymech.meshtools as mt

def add_h(s):
    return re.sub(r'(\d)(D\.rea)', r'_h\1\2', s)

def split_mesh(fname_in):
   meshfldr = '.'
   fname_out = add_h(fname_in)
   fname = os.path.join(meshfldr,fname_in)

   print('Read mesh '+fname)
   mesh = pm.readrea(fname)

   # cut mesh, only keep y > 0
   # also, get centroids for plot
   helems = []
   xyc  = np.empty((mesh.nel,2))
   for (iel, el) in enumerate(mesh.elem):
      x, y, _ = el.centroid
      xyc[iel,:] = x, y
      if x > 0.0:
         helems.append(iel)

   mt.keep_elements(mesh, helems, external_bc='SYM')

   xych = np.empty((mesh.nel,2))
   w = [ False for i in range(mesh.nel) ]
   s = [ False for i in range(mesh.nel) ]
   for (iel, el) in enumerate(mesh.elem):
      x, y, _ = el.centroid
      xych[iel,:] = x, y

   for (iel, el) in enumerate(mesh.elem):
      x, y, z = el.centroid
      xych[iel,:] = x, y
      for iface in range(2*mesh.ndim):
         bc = el.bcs[0, iface][0]
         fc = el.face_center(iface)
         #print(f'el={iel+1}, {iface}: bc = {bc}')
         if abs(fc[0]) < 1e-6:
            s[iel] = True
            el.bcs[0, iface][0] = 'SYM'
         if bc == 'MS':
            w[iel] = True
            el.bcs[0, iface][0] = 'W  '
         bc = el.bcs[0, iface][0]
         if bc == 'SYM' or bc == 'W  ':
            print(f'el={iel+1}, {iface}: bc = {bc}')
         

   fname = os.path.join(meshfldr,fname_out)
   print('Write mesh '+fname)
   pm.writerea(fname, mesh)
   print(f'new mesh={mesh.nel}')

   fig, ax = plt.subplots(1,1, figsize=(10, 8))
   ax.scatter(xyc[:,0], xyc[:,1], c = 'b', marker='o')
   ax.scatter(xych[:,0], xych[:,1], c = 'r', marker='x')
   ax.scatter(xych[w,0], xych[w,1], c = 'k', s = 150, marker='o')
   ax.scatter(xych[s,0], xych[s,1], c = 'g', s = 50, marker='o')

   plt.show()

if __name__ == "__main__":
   split_mesh('torus_coarse2D.rea')