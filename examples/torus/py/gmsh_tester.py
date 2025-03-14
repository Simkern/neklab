import sys, os, shutil, re
import subprocess

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

def prepare_test(mesh_params, fldr, basename_mesh, is_half):

    test_fldr = 'mesh_test'

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
            ffile=os.path.join(fldr,basename_mesh)
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
            print(f"\tSIZE file updated with lelg = {nel}.")

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
            print(f"\tSIZE file updated with nelf = {nelf}.")

            # Run mnl in the shell
            print("\trun makeneklab ...")
            # Open the log file for writing (this will overwrite the file)
            logfile = "build.txt"
            with open(logfile, "w") as log_file:
                # Run the command and pipe stdout and stderr to the logfile
                subprocess.run('makeneklab', stdout=log_file, stderr=log_file, cwd=tfldr, check=True)
            check_compilation(logfile)
        except subprocess.CalledProcessError as e:
            print(f"Error running mesh test: {e}")
            sys.exit(1)
        except FileNotFoundError as e:
            print(f"File not found error: {e}")
            sys.exit(1)
        
def test_mesh(fldr):
    test_fldr = 'mesh_test'
    # Check if the mesh_test folder exists
    tfldr = os.path.join(fldr, test_fldr)
    
    try:
        # Run nekbmpi with 'torus 12'
        print(f"\tRunning 'nekmpi torus 12'...")
        subprocess.run(['bash', '-c', 'nekmpi torus 12 > logfile.txt'], cwd=tfldr, check=True, stderr=subprocess.DEVNULL)
        pattern = r'Step'

        # Open the log file and read through each line
        with open(os.path.join(tfldr,'logfile.txt'), 'r') as file:
            for line in file:
                if re.search(pattern, line):
                    print(f"  {line.strip()}")
    
    except subprocess.CalledProcessError as e:
        print(f"\tError running 'nekmpi torus 12': {e}")
        sys.exit(1)