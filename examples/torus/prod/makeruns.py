import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import pandas as pd
import subprocess
import shutil
import argparse

def create_folder(home, fldr):
    if not os.path.exists(fldr):
        os.makedirs(fldr)
        print(f"Folder {fldr.replace(home+'/','')} created.")
        return True
    else:
        print(f"Folder {fldr.replace(home+'/','')} already exists. Skip.")
        return False

def copy_file(home, ifile, from_fldr, to_fldr):
   file = os.path.join(from_fldr,ifile)
   if os.path.exists(file):
      shutil.copy(file, to_fldr)
      print(f"\tcp: '{ifile}' -> {to_fldr.replace(home+'/','')}")
   else:
      print(f"Error: '{ifile}' not found. Skipping copy.")

def link_IC(home, from_fldr, sfile, to_fldr, dfile):
   src = os.path.join(from_fldr, sfile)
   dest = os.path.join(to_fldr, dfile)
   if os.path.exists(src):
      if not os.path.exists(dest):
         os.symlink(src, dest)
         print(f"\tslink: {sfile} -> {dest.replace(home+'/','')}")
      else:
         print(f"Symbolic link for {dfile} already exists in {to_fldr.replace(home+'/','')}. Skip.")
   else:
      print(f"Error: '{sfile}' not found in {from_fldr.replace(home+'/','')}. Skipping link creation.")

def link_mesh(home, geom_fldr, to_fldr, meshname):
   for ext in ['.re2', '.ma2']:
      sfile = meshname+ext
      dfile = 'torus'+ext
      src = os.path.join(geom_fldr, sfile)
      dest = os.path.join(to_fldr, dfile)
      if os.path.exists(src):
         if not os.path.exists(dest):
            os.symlink(src, dest)
            print(f"\tslink: {sfile} -> {dest.replace(home+'/','')}")
         else:
            print(f"Symbolic link for {sfile} already exists in {to_fldr.replace(home+'/','')}. Skip.")
      else:
         print(f"Error: '{sfile}' not found in {geom_fldr.replace(home+'/','')}. Skipping link creation.")

def check_compilation(logfile):
    success_msg = "#############################################################\n#                  Compilation successful!                  #\n#############################################################"
    # Open the log file in read mode
    try:
        with open(logfile, 'r') as f:
            content = f.read()
            if success_msg in content:
                print(success_msg)  # Print the success message if found
            else:
                print(f"Error in compilation. Abort.")
                sys.exit(1)
    except FileNotFoundError:
        print(f"Error: The file {logfile} was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(-1)

def compile_nek(compile_fldr):
   # Run mnl in the shell
   print("\tcompile nek ...")
   # Open the log file for writing (this will overwrite the file)
   logfile = "build.txt"
   with open(logfile, "w") as log_file:
      # Run the command and pipe stdout and stderr to the logfile
      subprocess.run('makeneklab', stdout=log_file, stderr=log_file, cwd=compile_fldr, check=True)
   check_compilation(logfile)
   print("\tdone.")

def write_usr_01(home, usr, usr_ref, data, tol, Tend):
   print(f'\nReading from\t{usr_ref.replace(home+'/','')}')
   print(f'Output to     \t{usr.replace(home+'/','')}\n')
   with open(usr_ref, 'r') as f:
      lines = f.readlines()
      with open(usr, 'w') as u:
         for line in lines:
            if 'real(dp), parameter :: womersley' in line:
               line = f'      real(dp), parameter :: womersley   = {data['Wo'].iloc[0]:.1f}_dp\n'
               print(line.strip())
            elif 'real(dp), parameter :: dp0ds' in line:
               line = f'      real(dp), parameter :: dp0ds       = {data['Re(fsamp_0)'].iloc[0]:14.12e}_dp\n'
               print(line.strip())
            elif 'real(dp), parameter :: dp1dsr' in line:
               line = f'      real(dp), parameter :: dp1dsr      = {data['Re(fsamp_1)'].iloc[0]:14.12e}_dp\n'
               print(line.strip())
            elif 'real(dp), parameter :: dp1dsi' in line:
               line = f'      real(dp), parameter :: dp1dsi      = {data['Im(fsamp_1)'].iloc[0]:14.12e}_dp\n'
               print(line.strip())
            elif 'real(dp), parameter :: tol' in line:
               line = f'      real(dp), parameter :: tol         = {tol:8.2e}_dp\n'
               print(line.strip())
            elif 'real(dp), parameter :: Tend' in line:
               line = f'      real(dp), parameter :: Tend        = {Tend:.1f}_dp\n'
               print(line.strip())
            u.write(line)
   return

def write_usr_02(home, usr, usr_ref, data, tol_nwt, tol_mf, tol_mode, max_nwt_iter):
   print(f'\nReading from\t{usr_ref.replace(home+'/','')}')
   print(f'Output to     \t{usr.replace(home+'/','')}\n')
   with open(usr_ref, 'r') as f:
      lines = f.readlines()
      with open(usr, 'w') as u:
         for line in lines:
            if 'real(dp), parameter :: womersley' in line:
               line = f'      real(dp), parameter :: womersley    = {data['Wo'].iloc[0]:.1f}_dp\n'
               print(line.strip())
            elif 'real(dp), parameter :: dp0ds' in line:
               line = f'      real(dp), parameter :: dp0ds        = {data['Re(fsamp_0)'].iloc[0]:14.12e}_dp\n'
               print(line.strip())
            elif 'real(dp), parameter :: dp1dsr' in line:
               line = f'      real(dp), parameter :: dp1dsr       = {data['Re(fsamp_1)'].iloc[0]:14.12e}_dp\n'
               print(line.strip())
            elif 'real(dp), parameter :: dp1dsi' in line:
               line = f'      real(dp), parameter :: dp1dsi       = {data['Im(fsamp_1)'].iloc[0]:14.12e}_dp\n'
               print(line.strip())
            elif 'real(dp), dimension(2), parameter :: mflow_target' in line:
               line = f'      real(dp), dimension(2), parameter :: mflow_target = [ 1.0_dp, {data['Q'].iloc[0]:5.3f}_dp ]\n'
               print(line.strip())
            elif 'real(dp) ::            tol_nwt' in line:
               line = f'      real(dp) ::            tol_nwt      = {tol_nwt:8.2e}_dp\n'
               print(line.strip())
            elif 'real(dp) ::            tol_mf' in line:
               line = f'      real(dp) ::            tol_mf       = {tol_mf:8.2e}_dp\n'
               print(line.strip())
            elif 'integer, parameter ::  tol_mode' in line:
               line = f'      integer, parameter ::  tol_mode     = {tol_mode:d}\n'
               print(line.strip())
            elif 'integer, parameter ::  max_nwt_iter' in line:
               line = f'      integer, parameter ::  max_nwt_iter = {max_nwt_iter:d}\n'
               print(line.strip())
            u.write(line)
   return

if __name__ == '__main__':

   parser = argparse.ArgumentParser()
   # Add the new mandatory arguments
   parser.add_argument('--Wo', type=float, required=True, help='Womersley number')
   parser.add_argument( '--Q', type=float, required=True, help='Pulsation amplitude.')
   parser.add_argument('--ref', type=str, help='Reference case')
   # Choose if the case should only be updated
   parser.add_argument('--update', action='store_true', help='Update instead of create.')
   # Choose if the data should be plotted for the relevant Wo
   parser.add_argument('--plot', action='store_true', help='Update instead of create.')
   
   # Parse arguments
   args = parser.parse_args()

   if len(sys.argv) == 5:
        update = False
        plot = False
   elif not args.update and not args.ref and not args.plot:
        # If any other argument is provided, print an error message and exit
        print("Error: Invalid argument provided.")
        print("Usage: p makeruns.py --Wo Wo --Q Q [--ref refdir] [--update] [--plot]")
        sys.exit(1)

   home      = os.getcwd()
   # parameters
   paramfile = os.path.join(home, 'Re_1700_delta_0.300.txt')
   # run
   run_fldr  = os.path.join(home, 'run')
   # compile
   cdir      = os.path.join(home, 'compile')
   # usr
   usr_ref   = 'torus.usr_ref'
   usr       = 'torus.usr'
   # geometry
   geom      = os.path.join(home, 'geom/half_torus')
   meshname  = 'torus_v1_3z_3D'
   # baseflow initial condition
   bfdir     = os.path.join(home, 'bf2d')

   Wo = args.Wo
   Q  = args.Q
   ReR          = 1700
   tol          = 1e-6
   Tend         = 200.0
   tol_nwt      = 1e-8
   tol_mf       = 1e-6 
   tol_mode     = 2
   max_nwt_iter = 10
   df0 = 0.0794

   params = pd.read_csv(paramfile, sep='\\s+', header=0)
   print('')
   print(params)

   select = (params['Wo'] == Wo)
   data = params[select]

   # Get the forcing columns
   fcols = params.columns[-3:]

   if Q in data['Q']:
      df = data[data['Q'] == Q]
   else:
      df_data = {
         'ReR': ReR,
         'Wo': Wo,
         'Q': Q,
         'T': 2*np.pi*ReR/Wo**2,
         'omega': Wo**2/ReR
      }
      for col in fcols:
         spl = CubicSpline(data['Q'], data[col])
         df_data[col] = spl(Q)
      df = pd.DataFrame([df_data])

   if Q >= data['Q'].max() or Q <= data['Q'].min():
      print('WARNING: Requested amplitude is outside the data range.')

   if args.plot:

      fig, axs = plt.subplots(1, 3, figsize=(15, 5))

      # Get the colormap 'jet'
      cmap = plt.get_cmap('jet')

      # Normalize Q values to map them to the colormap
      norm = plt.Normalize(vmin=params['Wo'].min(), vmax=params['Wo'].max())

      # Iterate through the last three columns and plot
      for i, col in enumerate(fcols):
         # Group by the unique values in 'Q'
         for wo_value in params['Wo'].unique():
            # Filter data for the current value of 'Q'
            subset = params[params['Wo'] == wo_value]
            
            # Plot the data for the current column
            color = cmap(norm(wo_value))  # Get color from colormap
            if wo_value == Wo:
               alpha = 1
               lw = 2
            else:
               alpha = 0.5
               lw = 1
            if wo_value % 5 == 0:
               axs[i].plot(subset['Q'], subset[col], 'o-', color=color, linewidth=lw, alpha=alpha, label=f'Wo={wo_value}')
            else:
               axs[i].plot(subset['Q'], subset[col], color=color, linewidth=lw, alpha=alpha)
         
         # Set labels and title for each subplot
         axs[i].set_xlabel('Q')
         axs[i].set_ylabel(col)
         axs[i].set_title(f'{col}')
         axs[i].scatter(Q, df[col], s=50, color='black', marker='x')

      axs[-1].legend()

      # Adjust the layout for a clean view
      plt.tight_layout()

      df_plt = np.abs(data[['Q', 'Re(fsamp_0)', 'Re(fsamp_1)', 'Im(fsamp_1)']] - [ 0.0, df0, 0.0, 0.0 ])
      df_plt.set_index('Q', inplace=True)
      ax = df_plt.plot(style=['o-', 'o-', 'o-'])
      ax.set_yscale('log')
      ax.set_xscale('log')
      for col in fcols:
         if 'amp_0' in col:
            ax.scatter(Q, df[col]-df0, s=50, color='black', marker='x')
         else:
            ax.scatter(Q, df[col], s=50, color='black', marker='x')
      plt.title(f'Wo = ${Wo:5.1f}$', fontsize=20)
      plt.xlabel('Q')
      plt.ylabel('f_0')
      plt.show()
      sys.exit()

   # Summary:
   print('Create/update the following folders:')
   dirnames = {}
   wdir = f'Wo_{Wo:05.1f}'
   dname = os.path.join(run_fldr, wdir)
   existdir = 'create.'
   if os.path.isdir(dname):
      existdir = f'exists.'
   msg = f'    {wdir}:'
   print(f'{msg:40s}{existdir:30s}')
   qdir = f'Q_{Q:5.3f}'
   casename = os.path.join(dname, qdir)
   existdir = 'create.'
   if os.path.isdir(casename):
      existdir = f'exists.'
   if not args.ref:
      msg = f'        {qdir}: startup'
   else:
      msg = f'        {qdir}: (ref:{args.ref})'
   print(f'{msg:40s}{existdir:30s}')
   refdir = None
   if args.ref:
      refdir = os.path.join(dname,args.ref)
      if not os.path.isdir(refdir):
         print(f'WARNING: reference dir {refdir} does not exist.')
      else:
         print(f'Refdir: {refdir}')

   print('\nSelected data:')
   print(df)

   if refdir is None:  # startup
      cpldir = os.path.join(cdir, '01_cold_start')
      ofile = os.path.join(cpldir,usr)
      rfile = os.path.join(cpldir,usr_ref)
      write_usr_01(home, ofile, rfile, df, tol, Tend)
   else: # newton
      cpldir = os.path.join(cdir, '02_newton')
      ofile = os.path.join(cpldir,usr)
      rfile = os.path.join(cpldir,usr_ref)
      write_usr_02(home, ofile, rfile, df, tol_nwt, tol_mf, tol_mode, max_nwt_iter)

   user_input = input(f"\nConfirm? (y/n) ")
   if user_input.lower() not in ['', 'y', 'yes']:
      sys.exit()

   compile_nek(cpldir)

   # create folders 
   create_folder(home, dname) # main folder
   if create_folder(home, casename) or args.update: # Q folders
      
      copy_file(home, 'nek5000', cpldir, casename)
      copy_file(home, 'torus.usr', cpldir, casename)
      
      if not args.update:
         copy_file(home, 'torus.par', cpldir, casename)
         link_mesh(home, geom, casename, meshname)
         if refdir is None:
            link_IC(home, bfdir, 'c2dtorus_v1_h_3z.fld', casename, 'c2d_rsttorus.fld')
         else:
            link_IC(home, refdir, 'BFNtorus0.f00001', casename, 'rsttorus0.f00001')