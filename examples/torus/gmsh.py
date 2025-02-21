import sys, argparse, os, json
import numpy as np
import matplotlib.pyplot as plt
from gmsh_writer import generate_mesh
from gmsh_tester import prepare_test, test_mesh
from gmsh_plotter import plot_gmsh
from read_2d_data import read_fields
from plot_2d_data import plot_2d_fld
from manipulate_2d_data import symmetrize_fld

geom_params = {
    'R': 1.0,
    'rt': 0.7,
    'rb': 0.7,
    'RBt': 0.95,
    'RBb': 0.95,
    'tht': np.pi/4.0,
    'thb': np.pi/4.0,
    'lambda1t': 0.8,
    'lambda1b': 0.8,
    'lambda2': 0.8,
    'dyc': 0.00
}

mesh_params = {
    'Lz': 1, 
    'Ncv': 11,
    'Nch': 11,
    'NB': 1,
    'NM': 5,
    'compressRatio_B': 0.9,
    'compressRatio_M': 0.87,
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
    fldr_h = 'geomh'
    param_fldr = 'mesh_params'
    basename2 = 'test2D'
    basename3 = 'test3D'

    params = {
        'geom_params': geom_params,
        'mesh_params': mesh_params
    }

    if args.plt:
        fig, ax = plt.subplots(1, 2, figsize=(20, 20))
        plot_gmsh(ax[0], geom_params, half=True, aux1=True, aux2=True, hlines=True)
        plot_gmsh(ax[1], geom_params, half=False, aux1=True, aux2=True, hlines=True)
        plt.show()
        sys.exit()
    else:
        if args.chk:
            fig, ax = plt.subplots(figsize=(20, 20))
            plot_gmsh(ax, geom_params, half=is_half, aux1=True, aux2=True, hlines=True)
            plt.show()

            user_input = input(f"Do you want to generate the GMSH script'? (y/n): ")
            is_generated_f, is_generated_h = False, False
            if user_input.lower() in ['', 'y', 'yes']:
                is_generated_f = generate_mesh(geom_params, mesh_params, fldr_h, basename2, is_half=True, confirm=True)
                is_generated_h = generate_mesh(geom_params, mesh_params, fldr,   basename2, is_half=False, confirm=True)
                with open(os.path.join(param_fldr,basename2+'.json'), 'w') as file:
                    json.dump(params, file, indent=4)
            if (is_generated_f, is_generated_h).any():

                testf = input("Generate mesh and test for full cross-section? [y]: ")
                if (testf.lower() in ['', 'y', 'yes'] and is_generated_f):
                    prepare_test(mesh_params, fldr,   basename3, is_half=False)
                    test_mesh(fldr)
                else:
                    print('Full mesh not tested.')
                
                testh = input("Generate mesh and test for half cross-section? [y]: ")
                if (testh.lower() in ['', 'y', 'yes'] and is_generated_f):
                    prepare_test(mesh_params, fldr_h, basename3, is_half=True)
                    test_mesh(fldr_h)
                else:
                    print('Half mesh not tested.')

        else:
            print('\nGenerate mesh and test for full cross-section:\n')
            if generate_mesh(geom_params, mesh_params, fldr,   basename2, is_half=False, confirm=False):
                prepare_test(mesh_params, fldr,   basename3, is_half=False)
                test_mesh(fldr)
            else:
                print('Full mesh not generated or tested.')
            print('\nGenerate mesh and test for half cross-section:\n')
            if generate_mesh(geom_params, mesh_params, fldr_h, basename2, is_half=True, confirm=False):
                prepare_test(mesh_params, fldr_h, basename3, is_half=True)
                test_mesh(fldr_h)
            else:
                print('Half mesh not generated or tested.')

    pattern = 'c2dtorus'
    runfldr   = os.path.join(fldr,   'mesh_test')
    runfldr_h = os.path.join(fldr_h, 'mesh_test')
    fname = pattern+'001.fld'
    if os.path.exists(os.path.join(runfldr,fname)) and os.path.exists(os.path.join(runfldr_h,fname)):
        # Create the figure and axis
        fig, ax = plt.subplots(1, 2, figsize=(20,8))
        #fig, ax2 = plt.subplots(1, 2, figsize=(20,8))
        ax[0].set_title('Full mesh')
        data, meta, nsteps = read_fields(pattern, cwd=runfldr)
        plot_2d_fld(ax[0], data.x, data.y, data.vx, draw_elements=True)
        #vxfs, vxfa = symmetrize_fld(x, y, vx)
        #plot_2d_fld(ax2[0], x, y, vxfs)
        #plot_2d_fld(ax2[1], x, y, vxfa)
        ax[1].set_title('Half mesh')
        data, meta, nsteps = read_fields(pattern, cwd=runfldr_h)
        plot_2d_fld(ax[1], data.x, data.y, data.vx, draw_elements=True)
        
        plt.show()
    else:
        print('Files not found')
        sys.exit()
