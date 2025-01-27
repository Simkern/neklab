import sys, argparse, os
import numpy as np
import matplotlib.pyplot as plt
from gmsh_writer import generate_mesh
from gmsh_tester import prepare_test, test_mesh
from gmsh_plotter import plot_gmsh
from read_2d_data import read_fields
from plot_2d_data import plot_2d_fld
from sort_2d_data import symmetrize_fld

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

    parser = argparse.ArgumentParser()
    parser.add_argument('--chk', action='store_true', help='Request user confirmation for each step.')
    parser.add_argument('--plt', action='store_true', help='Only plot the mesh structure.')

    # Parse arguments
    args = parser.parse_args()

    get_confirmation = False
    plot = False

    if len(sys.argv) == 1:
        get_confirmation = False
        plot = False
    elif not args.plt and not args.chk:
        # If any other argument is provided, print an error message and exit
        print("Error: Invalid argument provided.")
        print("Usage: p gmsh.py [chk|plt]")
        sys.exit(1)

    fldr   = 'geom'
    fldr_h = 'geom_h'
    basename2 = 'test2D'
    basename3 = 'test3D'

    if args.plt:
        plot_gmsh(geom_params, half=is_half, aux1=True, aux2=True, hlines=True)
    else:
        if args.chk:
            plot_gmsh(geom_params, half=is_half, aux1=True, aux2=True, hlines=True)

            user_input = input(f"Do you want to generate the GMSH script'? (y/n): ")
            if user_input.lower() in ['', 'y', 'yes']:
                generate_mesh(geom_params, mesh_params, fldr_h, basename2, is_half=True, confirm=True)
                generate_mesh(geom_params, mesh_params, fldr,   basename2, is_half=False, confirm=True)
            
            test_mesh_input = input(f"Do you want to test the mesh by running the necessary operations? (yes/no) [y]: ")
            if test_mesh_input.lower() in ['', 'y', 'yes']:
                prepare_test(mesh_params, fldr_h, basename3, is_half=True)
                prepare_test(mesh_params, fldr,   basename3, is_half=False)

            nekbmpi_input = input(f"Do you want to run the mesh tests? (yes/no) [y]: ")
            if nekbmpi_input.lower() in ['', 'y', 'yes']:
                test_mesh(fldr_h)
                test_mesh(fldr)
        else:
            print('\nGenerate mesh and test for full cross-section:\n')
            generate_mesh(geom_params, mesh_params, fldr,   basename2, is_half=False, confirm=False)
            prepare_test(mesh_params, fldr,   basename3, is_half=False)
            test_mesh(fldr)
            print('\nGenerate mesh and test for half cross-section:\n')
            generate_mesh(geom_params, mesh_params, fldr_h, basename2, is_half=True, confirm=False)
            prepare_test(mesh_params, fldr_h, basename3, is_half=True)
            test_mesh(fldr_h)

    pattern = 'c2dtorus'
    runfldr   = os.path.join(fldr,   'mesh_test')
    runfldr_h = os.path.join(fldr_h, 'mesh_test')
    fname = pattern+'001.fld'
    if os.path.exists(os.path.join(runfldr,fname)) and os.path.exists(os.path.join(runfldr_h,fname)):
        # Create the figure and axis
        fig, ax = plt.subplots(1, 2, figsize=(20,8))
        #fig, ax2 = plt.subplots(1, 2, figsize=(20,8))
        ax[0].set_title('Full mesh')
        x, y, vx, vy, vz, elmap, dt2d, metadata, nsteps = read_fields(pattern, cwd=runfldr)
        plot_2d_fld(ax[0], x, y, vx)
        #vxfs, vxfa = symmetrize_fld(x, y, vx)
        #plot_2d_fld(ax2[0], x, y, vxfs)
        #plot_2d_fld(ax2[1], x, y, vxfa)
        ax[1].set_title('Half mesh')
        x, y, vx, vy, vz, elmap, dt2d, metadata, nsteps = read_fields(pattern, cwd=runfldr_h)
        plot_2d_fld(ax[1], x, y, vx)
        
        plt.show()
    else:
        print('Files not found')
        sys.exit()