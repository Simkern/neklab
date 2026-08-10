import copy as cp
import numpy as np
import pymech as pm
import matplotlib.pyplot as plt

IC3D = 'IC3D_poiseuille.fld'
print(f'Read {IC3D}')
data3d = pm.neksuite.readnek(IC3D)

M2D = 'poiseuille.re2'
print(f'Read {M2D}')
m2d = pm.neksuite.readre2(M2D)

elem2d = []
elmap =  []
ctr    = []
for iel, elem in enumerate(data3d.elem):
   minz = elem.pos[2].min()
   if minz < 1e-6:
      elem2d.append(elem)
      elmap.append(iel)
      ctr.append(elem.centroid)

ctr3d = np.array(ctr)

data2d = cp.deepcopy(m2d)
data2d.var = (3, 3, 1, 0, 0)
data2d.time = 0.0
data2d.istep = 1
data2d.lr1 = (8, 8, 1)
data2d.elmap = np.array(elmap, dtype='int32')
for iel, elem in enumerate(elem2d):
   data2d.elem[iel].pos  = np.empty((3,8,8,1))
   data2d.elem[iel].vel  = np.empty((3,8,8,1))
   data2d.elem[iel].pres = np.empty((1,8,8,1))
   data2d.elem[iel].pos[:2,:,:,0] = elem.pos[:2,:,:,0]
   data2d.elem[iel].vel[:2,:,:,0] = elem.vel[:2,:,:,0]
   data2d.elem[iel].pres[0,:,:,0] = elem.pres[0,:,:,0]
   data2d.elem[iel].bcs = elem.bcs


pm.neksuite.writenek('IC2D_poiseuille.fld', data2d)
   

# extract first slice