import os, sys, re
import shutil
import subprocess
import numpy as np

def update_viscosity(filename, Re_value):
   # Read the input file
   with open(filename, 'r') as file:
      file_contents = file.read()
   # Replace the viscosity line with the new value corresponding to -Re
   file_contents = re.sub(r'viscosity = -40\.0', f'viscosity = {-Re_value}', file_contents)
   # Save the modified contents to the output file
   with open(filename, 'w') as file:
      file.write(file_contents)

def create_runfolder(step, base):
   # Create a new folder inside the newton folder
   new_folder= os.path.join(base, run_folder_base + f'{step:02d}')
   if not os.path.exists(new_folder):
      os.makedirs(new_folder)
      print(f"\tCreated directory: {new_folder}")
      skip = False
   else:
      print(f"\tDirectory exists: {new_folder}")
      skip = True
      return new_folder, skip
 
   # Copy the input par file to the new folder
   par_file = os.path.join(base, '1cyl.par')
   shutil.copy(par_file, new_folder)
   print(f"\tCopied {par_file} to {new_folder}")
   shutil.copy('run_nek.sh', new_folder)
 
   # Create symbolic links for the .re2 and nek5000 files in the new folder
   link_files = [ '1cyl.re2', '1cyl.ma2', 'nek5000' ]
   for file in link_files:
      os.symlink(os.path.join('..', file), os.path.join(new_folder, file))
   print(f"\tCreated symbolic links in {new_folder}")

   return new_folder, skip

def copy_from_to(from_file, to_file):
   skip = False
   if os.path.exists(from_file):
      if os.path.exists(to_file):
         #print(f"{to_file} exists. Do nothing.")
         skip = True
      else:
         shutil.copy(from_file, to_file)
         print(f"\tCopied {from_file} to {to_file}")
   else:
      print(f"\tError: File {from_file} not found.")
      sys.exit()
   return skip

def launch_nek():
   try:
      print("Running simulation...")
      subprocess.run(['sh', 'run_nek.sh'], check=True, shell=False)
      print("Simulation finished successfully.")
   except subprocess.CalledProcessError as e:
      print(f"Simulation failed with error: {e}")
      sys.exit()
   return

def compute_iteration(Re, step):
   # Create a new folder inside the newton folder
   new_folder, skip = create_runfolder(step, folder)
   
   if (not skip):
      # Copy the initial condition from the appropriate location
      if step == 1:
         from_file = os.path.join(folder, 'BF.fld')
         to_file = os.path.join(new_folder, 'BF.fld')
         skiprun = copy_from_to(from_file, to_file)
      else:
         from_file = bf_base_file = f'BF{step-1:02d}.fld'
         to_file = os.path.join(new_folder, 'BF.fld')
         skiprun = copy_from_to(from_file, to_file)
         from_file = bf_base_file = f'EV{step-1:02d}.fld'
         to_file = os.path.join(new_folder, 'EV.fld')
         skiprun = copy_from_to(from_file, to_file)

      if step > 1:
         par_file = os.path.join(new_folder, '1cyl.par')
         update_viscosity(par_file, Re)
      
      # Change the directory to the new run folder
      os.chdir(new_folder)
      # Execute the simulation command: nekmpi 1cyl 12
      launch_nek()
      # return home
      os.chdir(home)
 
   # Step 3: Copy the new results
   # baseflow
   from_file = os.path.join(new_folder, 'BF_1cyl0.f00001')
   to_file = f'BF{step:02d}.fld'
   skip = copy_from_to(from_file, to_file)
   # eigenvector
   from_file = os.path.join(new_folder, 'dir1cyl0.f00001')
   to_file = f'EV{step:02d}.fld'
   skip = copy_from_to(from_file, to_file)
   # eigenvalues
   from_file = os.path.join(new_folder, 'dir_eigenspectrum.npy')
   to_file = f'dir_eigenspectrum{step:02d}.npy'
   skip = copy_from_to(from_file, to_file)

   return to_file

if __name__ == "__main__":

   # Define paths
   newton_folder    = 'newton'
   stability_folder = 'stability'
   folder           = 'newton_eigs'
   run_folder_base = 'run_'
   par_file = '1cyl.par'

   unstable = False
   Re = 40.0
   Rev = [ Re ]
   Re_inc = 2.0
   step = 0

   home = os.getcwd()
   print(f'CWD: {home}\n')

   while (not unstable):
      step += 1

      if step > 1:
         Re  += Re_inc
         Rev.append(Re)

      eig_file = compute_iteration(Re, step)

      data = np.load(eig_file)
      print(Re)
      print(data)

      if data[:,0].max() > 0.0:
         unstable = True

Re_a, Re_b = Rev[-2:]
Re_d = Re_b - Re_a

print('\nStart bisection:\n')
# Golden ratio bisection
golden_ratio = (1 + 5 ** 0.5) / 2  # Golden ratio constant
tol = 1e-4
max_iter = 12

do_a, do_b = True, True

while Re_d > tol and step < max_iter:
   Re1 = Re_b - Re_d / golden_ratio
   Re2 = Re_a + Re_d / golden_ratio
   print('bisection', Re1, Re2)

   if do_a:
      step += 1
      print(f'\n\tNext Re = {Re1}\n')
     
      eig_file = compute_iteration(Re1, step)
      
      data1 = np.load(eig_file)
      print(data1)
      f1 = abs(data1[:,0].max())

   if do_b:
      step += 1
      print(f'\n\tNext Re = {Re2}\n')
      
      eig_file = compute_iteration(Re2, step)

      data2 = np.load(eig_file)
      print(data2)
      f2 = abs(data2[:,0].max())

   do_a, do_b = False, False
   # Narrow the search interval based on function values
   if f1 < f2:
       Re_b = Re2  # Minimum is in the left part, so adjust Re_b
       Rev.append(Re_b)
       do_a = True
   else:
       Re_a = Re1  # Minimum is in the right part, so adjust Re_a
       Rev.append(Re_a)
       do_b = True
   Re_d = Re_b - Re_a