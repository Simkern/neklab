import os, sys, glob, copy
import numpy as np
import matplotlib.pyplot as plt
import pymech as pm

from check_mesh import *

mdir = 'geom'

m2Dname = 'torus_v2_2D.re2'
m3Dname = 'torus_v2_3z_3D.re2'

m2Dfile = os.path.join(mdir,m2Dname)
m3Dfile = os.path.join(mdir,m3Dname)

m2Df = pm.neksuite.readre2(m2Dfile)
m3Df = pm.neksuite.readre2(m3Dfile)

m3Df_s = repair_extruded_mesh(m2Df, m3Df, plot=True)

m3Dname_out = os.path.join(mdir,m3Dname.replace('.re2','_c.re2'))

pm.neksuite.writere2(m3Dname_out, m3Df_s)

