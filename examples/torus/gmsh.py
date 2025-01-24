import sys, os, shutil, time
import subprocess
import numpy as np
from gmsh_writer import generate_gmsh_script
from gmsh_plotter import plot_gmsh

geom_params = {
    'R': 1.0,
    'rt': 0.5,
    'rb': 0.7,
    'RBt': 0.92,
    'RBb': 0.94,
    'tht': np.pi/6.0,
    'thb': np.pi/3.0,
    'lambda1t': 0.75,
    'lambda1b': 0.6,
    'lambda2': 0.5,
    'dyc': 0.04
}

mesh_params = {
    'Lz': 1, 
    'Ncv': 11,
    'Nch': 11,
    'NB': 1,
    'NM': 5,
    'compressRatio_B': 0.6,
    'compressRatio_M': 0.95,
    'Nz': 180
}

is_half = True

if __name__ == "__main__":

    test_fldr = 'mesh_test'
    if is_half:
        fldr = 'geom_h'
    else:
        fldr = 'geom'
    mkscript = 'mkmsh.sh'
    basename = 'test2D'
    geoname = basename+'.geo'
    filename = os.path.join(fldr,geoname)

    plot_gmsh(geom_params, half=is_half, aux1=True, aux2=True, hlines=True)
   
    user_input = input(f"Do you want to generate the GMSH script with the filename '{filename}'? (y/n): ")
    if user_input.lower() in ['', 'y', 'yes']:

        if os.path.exists(filename):
            overwrite_input = input(f"The file '{filename}' already exists. Do you want to overwrite it? (yes/no) [y]: ")
            if overwrite_input.lower() not in ['', 'y', 'yes']:
                print("File will not be overwritten. Aborting script generation.")
                sys.exit(0)

        generate_gmsh_script(geom_params, mesh_params, half=is_half, meshDim=2, filename=filename)

        # Prompt to run GMSH
        run_input = input(f"Do you want to visualize the mesh in '{filename}'? (y/n) [n]: ")
        
        if run_input.lower() in ['y', 'yes']:
            try:
                # Run GMSH with the generated script
                subprocess.run(['gmsh', filename], check=True)
                print(f"GMSH ran successfully with the script '{filename}'.")
            except subprocess.CalledProcessError as e:
                print(f"Error running GMSH: {e}")

        # Prompt to run GMSH
        run_input = input(f"Do you want to save the GMSH mesh from script '{filename}'? (y/n) [y]: ")
        
        if run_input.lower() in ['', 'y', 'yes']:
            try:
                # Run GMSH with the generated script
                subprocess.run(['gmsh', filename, '-2'], check=True)
                print(f"GMSH ran successfully with the script '{filename}'.")
            except subprocess.CalledProcessError as e:
                print(f"Error running GMSH: {e}")

        # Ask whether to prepare the mesh
        prepare_mesh_input = input(f"Do you want to prepare the mesh to run by executing 'bash {mkscript} {basename}'? (yes/no) [y]: ")
        if prepare_mesh_input.lower() in ['', 'y', 'yes']:
            # Check if the folder exists
            if os.path.exists(fldr):
                try:
                    # Change directory to folder and run the bash script
                    subprocess.run(['bash', mkscript, basename], cwd=fldr, check=True)
                    print(f"Mesh preparation script {mkscript} ran successfully with '{basename}'.")
                except subprocess.CalledProcessError as e:
                    print(f"Error running {mkscript}: {e}")
        
        test_mesh_input = input(f"Do you want to test the mesh by running the necessary operations? (yes/no) [y]: ")
        if test_mesh_input.lower() in ['', 'y', 'yes']:
            # Check if the mesh_test folder exists
            tfldr = os.path.join(fldr, test_fldr)
            Nh = mesh_params['Nch'] - 1
            Nv = mesh_params['Ncv'] - 1
            NB = mesh_params['NB']
            NM = mesh_params['NM'] - 1
            nelf = Nh*Nv + 2*NM*(Nh+Nv) + 2*NB*(Nh+Nv)
            if is_half:
                nelf = int(nelf/2.0)
            nel = 5*nelf
            if os.path.exists(tfldr):
                try:
                    # Navigate to mesh_test folder
                    # Verbose copy from ../test3D.re2 and ../test3D.ma2 to torus.re2 and torus.ma2
                    ffile=os.path.join(fldr,'test3D')
                    tfile=os.path.join(tfldr,'torus')
                    for ext in ['.re2', '.ma2']:
                        shutil.copy2(ffile+ext, tfile+ext)
                        print(f"{ffile+ext} -> {tfile+ext}")

                    # Modify the SIZE file to replace 'lelg=540' with new value
                    ffile=os.path.join(tfldr,'SIZE_ref')
                    tfile=os.path.join(tfldr,'SIZE')
                    with open(ffile, 'r') as file:
                        content = file.read()
                    content = content.replace('parameter (lelg=540)', f'parameter (lelg={nel})')
                    with open(tfile, 'w') as file:
                        file.write(content)
                    print(f"SIZE file updated with lelg = {nel}.")

                    # Modify the torus.usr file to replace 'nelf=54' with new value
                    ffile=os.path.join(tfldr,'torus.usr_ref')
                    tfile=os.path.join(tfldr,'torus.usr')
                    with open(ffile, 'r') as file:
                        content = file.read()
                    content = content.replace('integer, parameter :: nelf     = 54', f'integer, parameter :: nelf     = {nelf}')
                    if is_half:
                        content = content.replace("!call setbc(2,1,'SYM')", f"call setbc(2,1,'SYM')")
                        content = content.replace('logical, parameter :: if_sym   = .false.', f'logical, parameter :: if_sym   = .true.')
                    with open(tfile, 'w') as file:
                        file.write(content)
                    print(f"SIZE file updated with nelf = {nelf}.")

                    # Run mnl in the shell
                    print("run makeneklab ...")
                    subprocess.run(['makeneklab > build.txt'], cwd=tfldr, check=True)
                    print("done.")
                except subprocess.CalledProcessError as e:
                    print(f"Error running mesh test: {e}")
                except FileNotFoundError as e:
                    print(f"File not found error: {e}")

        nekbmpi_input = input(f"Do you want to run the test? (yes/no) [y]: ")
        if nekbmpi_input.lower() in ['', 'y', 'yes']:
            try:
                # Run nekbmpi with 'torus 12'
                subprocess.run(['bash', '-c', 'nekmpi torus 12'], cwd=tfldr, check=True)
                print("Running 'nekmpi torus 12'...")

            except subprocess.CalledProcessError as e:
                print(f"Error running 'nekmpi torus 12': {e}")
    