import os, sys
import subprocess
import shutil

# Define the tol range and tol_mode values
tol_values = [1e-12, 1e-10, 1e-8, 1e-6]
tol_mode_values = [1, 2]

# Path to the .usr file
usr_file_path = "1cyl.usr"

# Path to the base directory and BF.fld
base_dir = "01_run_base"
bf_fld_path = os.path.join(base_dir, "BF.fld")

# Directory for storing the run cases
run_base_dir = "run_case"

# Create the base directory for the run cases if it doesn't exist
if not os.path.exists(run_base_dir):
    os.makedirs(run_base_dir)

# Path to the utils directory where "mkrun" is located
utils_dir = "/home/skern/utils"

# Source the .bash_aliases to use the alias for mnlc
bashrc = os.path.expanduser("~/.bashrc")

# Iterate through the combinations of tol and tol_mode
case_counter = 1
for tol in tol_values:
    for tol_mode in tol_mode_values:
        # 1. Replace the corresponding lines in the .usr file
        with open(usr_file_path, 'r') as usr_file:
            usr_lines = usr_file.readlines()
        
        # Update the lines with tol and tol_mode (adjust the line indices based on your .usr file structure)
        for i, line in enumerate(usr_lines):
            if "tol = " in line:
                usr_lines[i] = f"         tol = {tol}\n"  # Update tol
            elif "tol_mode = " in line:
                usr_lines[i] = f"      tol_mode = {tol_mode}\n"  # Update tol_mode

        # Write the updated content back to the .usr file
        with open(usr_file_path, 'w') as usr_file:
            usr_file.writelines(usr_lines)

        # 2. Run "mnlc" (replace with the actual command to run your Fortran simulation)
        print(f"Running mnl for case {case_counter} with tol={tol} and tol_mode={tol_mode}")
        subprocess.run(f"bash -i -c 'source {bashrc} && snlt && mnl'", shell=True, check=True)

        # 3. Create the foldername and make the directory
        foldername = f"{case_counter:02d}_run_newton_test"
        print(foldername)

        if not os.path.exists(foldername):
            # 4. Run "mkrun 1cyl fldr"
            print(f"Running mkrun for folder {foldername}")
            subprocess.run(f"bash -i -c 'source {bashrc} && {utils_dir}/mkrun 1cyl {foldername}'", shell=True, check=True)
        else:
            # 4. Run "udrun 1cyl fldr"
            print(f"Running udrun for folder {foldername}")
            subprocess.run(f"bash -i -c 'source {bashrc} && {utils_dir}/udrun 1cyl {foldername}'", shell=True, check=True)

        summary_file_path = os.path.join(foldername, "setup_summary.txt")
        with open(summary_file_path, 'w') as summary_file:
            summary_file.write(f"Setup Summary for Case {case_counter}\n")
            summary_file.write("="*30 + "\n")
            summary_file.write(f"tol = {tol}\n")
            summary_file.write(f"tol_mode = {tol_mode}\n")
            summary_file.write(f"Input .usr file: {usr_file_path}\n")
            summary_file.write(f"Base directory: {base_dir}\n")
            summary_file.write(f"BF.fld path: {bf_fld_path}\n")
            summary_file.write("="*30 + "\n")
        
        # 5. Copy the BF.fld file into the new folder
        print(f"Copying BF.fld into {foldername}")
        shutil.copy(bf_fld_path, foldername)

        # Increment the case counter
        case_counter += 1

print("All cases completed!")
