import sys, os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('/ccc/work/cont003/gen6362/kernjoha/projects/neklab_torus/prod/py')
from manipulate_2d_data import h2d_to_f2d

def get_wdir(Wo, prefix=None):
   return f'{prefix}Wo_{Wo:5.1f}'

def get_qdir(Q):
   return f'Q_{Q:5.3f}'

if __name__ == "__main__":

   home = '/ccc/scratch/cont003/gen6362/kernjoha/projects/neklab_torus/prod'
   rundir = 'v0_stab'

   # mesh version 0
   f2d_file_ref = 'bf2d/c2dtorus_v0_f.fld'
   h2d_pattern  = 'n2dtorus'

   Wo = 35.
   Q = 0.001

   hfldr   = os.path.join(home, get_wdir(Wo), get_qdir(Q))
   outfldr = os.path.join(home, rundir)

   # Create the full version of half the torus for a series of flds
   h2d_to_f2d(h2d_pattern=h2d_pattern, f2d_file_ref=f2d_file_ref, hfldr=hfldr, outfldr=outfldr, outpattern='f2dtorus')