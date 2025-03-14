import numpy as np
import matplotlib.pyplot as plt

def compute_bisector(xt, yt, xb, yb):
   # connector
   dx =  xb -xt
   dy = -yb -yt
   norm = np.sqrt(dx**2 + dy**2)
   # bisector
   ux  =  dy/norm
   uy  = -dx/norm
   ul = np.array([-ux, uy])
   ur = np.array([ ux, uy])
   return ul, ur, norm

def compute_auxiliary_geometry_data(rings):

    pts_aux   = []
    midpoints = []

    for ring in rings:
        xt, yt, xb, yb, length = ring

        ## INNER RING
        # midpoint
        xm = (xt + xb)/2.0
        ym = (yt - yb)/2.0
        # bisector
        _, ur, norm = compute_bisector(xt, yt, xb, yb)
        ux  = ur[0]            # x-component of normalized bisector
        uy  = ur[1]            # y-component of normalized bisector
        h   = norm/2.0         # half-distance between arc points
        L   = np.sqrt(length**2 - h**2) # distance along bisector
        # auxiliary points
        midpoints.append([ xm, ym ])
        midpoints.append([-xm, ym ])
        pts_aux.append([ -xm - L*ux, ym + L*uy ])
        pts_aux.append([  xm + L*ux, ym + L*uy ])

    return pts_aux, midpoints