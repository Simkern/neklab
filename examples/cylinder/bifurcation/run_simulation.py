import os, re
import shutil
import subprocess

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
      print(f"Created directory: {new_folder}")

   print(os.getcwd())
 
   # Copy the input par file to the new folder
   par_file = os.path.join(base, '1cyl.par')
   shutil.copy(par_file, new_folder)
   print(f"Copied {par_file} to {new_folder}")
   shutil.copy('../run_nek.sh', new_folder)
 
   # Create symbolic links for the .re2 and nek5000 files in the new folder
   link_files = [ '1cyl.re2', '1cyl.ma2', 'nek5000' ]
   for file in link_files:
      os.symlink(os.path.join('..', file), os.path.join(new_folder, file))
   print(f"Created symbolic links in {new_folder}")

   return new_folder

if __name__ == "__main__":

   # Define paths
   newton_folder    = 'newton'
   stability_folder = 'stability'
   run_folder_base = 'run_'
   par_file = '1cyl.par'

   unstable = False
   Re = 40.0
   Re_inc = 2.0
   step = 0

   home = os.getcwd()

   while (not unstable and Re < 41.):
      step += 1

      if step > 1:
         Re += Re_inc

      # Create a new folder inside the newton folder
      new_folder_newton = create_runfolder(step, newton_folder)

      if step > 1:
         par_file = os.path.join(new_folder_newton, '1cyl.par')
         update_viscosity(par_file, Re)

      # Change the directory to the new run folder
      os.chdir(new_folder_newton)

      # Execute the simulation command: nekmpi 1cyl 12
      try:
         print("Running simulation...")
         subprocess.run(['sh', 'runnek.sh'], check=True, shell=False)
         print("Simulation finished successfully.")
      except subprocess.CalledProcessError as e:
         print(f"Simulation failed with error: {e}")

      os.chdir(home)

      # Step 2: After simulation completes, create the same subfolder inside the stability folder
      new_folder_stability = create_runfolder(step, stability_folder)

      if step > 1:
         par_file = os.path.join(new_folder_stability, '1cyl.par')
         update_viscosity(par_file, Re)

      # Step 3: Copy the result file from the newton folder to the stability folder
      nwt_output_file = os.path.join(new_folder_newton, 'nwt1cyl0.f00001')
      bf_fld_file     = os.path.join(new_folder_stability, 'BF.fld')

      if os.path.exists(nwt_output_file):
         shutil.copy(nwt_output_file, bf_fld_file)
         print(f"Copied {nwt_output_file} to {bf_fld_file}")
      else:
         print(f"Error: File {nwt_output_file} not found. Ensure the simulation has run successfully.")

      # Execute the simulation command: nekmpi 1cyl 12
      #try:
      #   print("Running simulation...")
      #   subprocess.run(["nekmpi", "1cyl", "12"], check=True)
      #   print("Simulation finished successfully.")
      #except subprocess.CalledProcessError as e:
      #   print(f"Simulation failed with error: {e}")