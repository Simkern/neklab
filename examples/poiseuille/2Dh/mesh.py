import copy as cp
import numpy as np
import pymech as pm
import matplotlib.pyplot as plt

M2D = 'IC_poiseuille_uv_r.fld'
print(f'Read {M2D}')
m2d = pm.neksuite.readnek(M2D)

x2d = np.array([ e.pos[0].ravel() for e in m2d.elem ]).ravel()
y2d = np.array([ e.pos[1].ravel() for e in m2d.elem ]).ravel()
xy  = np.stack((x2d, y2d)).T
#xy  = np.unique(np.stack((x2d, y2d)).T, axis=0)
c   = np.array([ e.centroid for e in m2d.elem ] )

fig, ax = plt.subplots(1,1, figsize=(12, 10))
z0 = 0.0
z1 = np.pi/2
ax.scatter(x2d, y2d)
ax.scatter(c[:,0], c[:,1])

with open('poiseuille.his', 'w') as f:
   f.write(f'{2*xy.shape[0]:d}\n')
   for (x,y) in xy:
      f.write(f'{x:.8f} {y:.8f} {z0:.8f}\n')
   for (x,y) in xy:
      f.write(f'{x:.8f} {y:.8f} {z1:.8f}\n')

plt.show()