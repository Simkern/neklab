import sys, os, glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import pandas as pd
import subprocess
import shutil
import argparse

if __name__ == '__main__':

   focus = 25
   
   paramfile = 'Re_1700_delta_0.300.txt'
   params = pd.read_csv(paramfile, sep='\\s+', header=0)
   print(params[params['Wo']==100].to_string(float_format='%.12f'))

   sys.exit()

   dbfile = 'db_Nek5000.csv'
   df = pd.read_csv(dbfile)
   df.set_index(['Wo', 'Q'], inplace=True)

   db10file = 'db_Nek5000_P10.csv'
   df10 = pd.read_csv(db10file)
   df10.set_index(['Wo', 'Q'], inplace=True)

   fig, axs = plt.subplots(1, 3, figsize=(15, 5))

   cnames = [ 'Re(dpds_00)', 'Re(dpds_01)', 'Im(dpds_01)' ]
   fcols = params.columns[-3:]
   # Get the colormap 'jet'
   cmap = plt.get_cmap('jet')

   # Normalize Q values to map them to the colormap
   Wov = np.unique(df.index.get_level_values('Wo'))
   Wov10 = np.unique(df10.index.get_level_values('Wo'))
   norm = plt.Normalize(vmin=Wov.min(), vmax=Wov.max())

   # Iterate through the last three columns and plot
   for wo_value in Wov:
      # Filter data for the current value of 'Q'
      # 
      subset = df.xs(wo_value, level='Wo')
      if wo_value in Wov10:
         sub10  = df10.xs(wo_value, level='Wo')
      subp   = params[params['Wo'] == wo_value]
      
      # Plot the data for the current column
      color = cmap(norm(wo_value))  # Get color from colormap
      alpha = 1
      if wo_value == focus:
          lw = 3
      else:
          lw = 1
      
      for i, (n,f) in enumerate(zip(cnames,fcols)):
         if i == 0:
            sh1 = 0.07948553
            sh2 = 0.07945995
            if sh1 == 0.0:
                axs[i].plot(subset.index, subset[n], 'o-', color=color, linewidth=lw, alpha=alpha)
                if wo_value in Wov10:
                   axs[i].plot(sub10.index, sub10[n], ':', color=color, linewidth=lw*2, alpha=alpha)
            else:
                axs[i].plot(subset.index, (subset[n]-sh1)/sh1*100, 'o-', color=color, linewidth=lw, alpha=alpha)
                if wo_value in Wov10:
                   axs[i].plot(sub10.index, (sub10[n]-sh1)/sh1*100, ':', color=color, linewidth=lw*2, alpha=alpha)
            if sh2 == 0.0:
                axs[i].plot(subp['Q'], subp[f], '--', color=color, linewidth=lw*2, alpha=alpha)
            else:
                axs[i].plot(subp['Q'], (subp[f]-sh2)/sh2*100, '--', color=color, linewidth=lw*2, alpha=alpha)
            
         else:
            if i == 2 and wo_value % 5 == 0:
               axs[i].plot(subset.index, subset[n], 'o-', color=color, linewidth=lw, alpha=alpha, label=f'Wo={wo_value}')
               axs[i].plot(subp['Q'], subp[f], '--', color=color, linewidth=lw*2, alpha=alpha, label=f'Wo={wo_value}, M')
            else:
               axs[i].plot(subset.index, subset[n], 'o-', color=color, linewidth=lw, alpha=alpha)
               axs[i].plot(subp['Q'], subp[f], '--', color=color, linewidth=lw*2, alpha=alpha)
      
   # Set labels and title for each subplot
   
   for i, ax in enumerate(axs):
      ax.set_xlabel('Q')
      ax.set_ylabel(cnames[i])
      if i == 0:
         tle = f'({cnames[i]} - f_steady)/f_steady (%)'
      else:
         tle = f'{cnames[i]}'
      ax.set_title(tle)

   axs[-1].legend()

   qv = np.array([ 0.001, 0.002, 0.003, 0.005, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100 ])
   # Get the colormap 'jet'
   cmap = plt.get_cmap('jet')
   # Normalize Q values to map them to the colormap
   norm = plt.Normalize(vmin=qv.min(), vmax=qv.max())
   fig, axs = plt.subplots(1, 3, figsize=(15, 5))

   for qp in qv:
      color = cmap(norm(qv))  # Get color from colormap
      subset = df.xs(qp, level='Q')
      for i, n in enumerate(cnames):
         axs[i].plot(subset.index, subset[n], 'o-', label=f'Q={qp}')

   for i, ax in enumerate(axs):
      ax.set_xlabel('Wo')
      ax.set_ylabel(cnames[i])
      tle = f'{cnames[i]}'
      ax.set_title(tle)

   axs[-1].legend()

   # Adjust the layout for a clean view
   plt.tight_layout()


   plt.show()
