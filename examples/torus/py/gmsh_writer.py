import sys, os, shutil, time
import subprocess
import numpy as np
from gmsh_compute import compute_auxiliary_geometry_data

# Function to create Point
def create_point(ip, a, b, c, d):
    return f"Point({ip}) = {{{a}, {b}, {c}, {d}}};"

# Function to create Circle
def create_circle(icl, start, center, end):
    return f"Circle({icl})={{{start}, {center}, {end}}};"

# Function to create Line
def create_line(icl, start, end):
    return f"Line({icl})={{ {start}, {end} }};"

# Function to create Transfinite Line
def create_transfinite_line(lines, nc, progression=None, bump=None):
    if progression:
        return f"Transfinite Line {{ {', '.join(map(str, lines))} }} = {nc} Using Progression {progression};"
    elif bump:
        return f"Transfinite Line {{ {', '.join(map(str, lines))} }} = {nc} Using Bump {bump};"
    else:
        return f"Transfinite Line {{ {', '.join(map(str, lines))} }} = {nc};"

# Function to create Line Loop and Plane Surface
def create_line_loop_surface(il, lines, surface_name):
    return f"Line Loop({il})={{ {', '.join(map(str, lines))} }};   Plane Surface({surface_name})={{ {surface_name} }};"

def generate_mesh(geom_params, mesh_params, fldr, basename, is_half, confirm=False):

    mkscript = 'mkmsh.sh'
    geoname = basename+'.geo'
    filename = os.path.join(fldr,geoname)
    is_generated = False

    if os.path.exists(filename):
        overwrite_input = input(f"The file '{filename}' already exists. Do you want to overwrite it? (yes/no) [y]: ")
        if overwrite_input.lower() not in ['', 'y', 'yes']:
            print("\tFile will not be overwritten. Aborting script generation.")
            return is_generated

    generate_gmsh_script(geom_params, mesh_params, half=is_half, meshDim=2, filename=filename)
    is_generated = True

    if confirm:
        # Prompt to run GMSH
        run_input = input(f"Do you want to visualize the mesh in '{filename}'? (y/n) [n]: ")
    else:
        run_input = 'n'
    
    if run_input.lower() in ['y', 'yes']:
        try:
            # Run GMSH with the generated script
            subprocess.run(['gmsh', filename], check=True)
            print(f"\tGMSH ran successfully with the script '{filename}'.")
        except subprocess.CalledProcessError as e:
            print(f"Error running GMSH: {e}")
            sys.exit(1)

    if confirm:
        # Prompt to run GMSH
        run_input = input(f"Do you want to save the GMSH mesh from script '{filename}'? (y/n) [y]: ")
    else:
        run_input = 'y'
    
    if run_input.lower() in ['', 'y', 'yes']:
        try:
            # Run GMSH with the generated script
            with open("out_gmsh.txt", "w") as log_file:
                # Run the command and pipe stdout and stderr to the logfile
                subprocess.run(['gmsh', filename, '-2 -order 2'], stdout=log_file, stderr=log_file, check=True)
            print(f"\tGMSH ran successfully with the script '{filename}'.")
        except subprocess.CalledProcessError as e:
            print(f"Error running GMSH: {e}")
            sys.exit(1)

    if confirm:
        # Ask whether to prepare the mesh
        prepare_mesh_input = input(f"Do you want to prepare the mesh to run by executing 'bash {mkscript} {basename}'? (yes/no) [y]: ")
    else:
        prepare_mesh_input = 'y'

    if prepare_mesh_input.lower() in ['', 'y', 'yes']:
        # Check if the folder exists
        if os.path.exists(fldr):
            try:
                # Change directory to folder and run the bash script
                subprocess.run(['bash', mkscript, basename], cwd=fldr, check=True)
                print(f"\tMesh preparation script {mkscript} ran successfully with '{basename}'.")
            except subprocess.CalledProcessError as e:
                print(f"Error running {mkscript}: {e}")
                sys.exit(1)

    return is_generated

def generate_gmsh_script(geom, mesh, half=False, meshDim=2, filename="pipe_mesh.geo"):
    # Extract data
    R, rt, rb, RBt, RBb, tht, thb, lambda1t, lambda1b, lambda2, dyc = geom.values()
    Lz, Nch, Ncv, NB, NM, compressRatio_B, compressRatio_M, Nz, = mesh.values()

    # compute data
    ra = (rt + rb)/2.0
    #top
    cost = np.cos(tht)
    sint = np.sin(tht)
    dxt  = rt * cost
    dyt  = rt * sint
    dxBt = RBt * cost
    dyBt = RBt * sint
    # bottom
    cosb = np.cos(thb)
    sinb = np.sin(thb)
    dxb  = rb * cosb
    dyb  = rb * sinb
    dxBb = RBb * cosb
    dyBb = RBb * sinb

    rings = []
    rings.append([ dxt,  dyt,  dxb,  dyb, lambda2*R + ra])
    rings.append([dxBt, dyBt, dxBb, dyBb, RBb])
    
    pts_aux, _ = compute_auxiliary_geometry_data(rings)

    # General Settings Section
    general_settings = f"""
/*
   ***** gmsh Script for generating 2D/3D mesh for straight pipe. *****
   ** Saleh Rezaeiravesh, salehr@kth.se
   >>> For nomenclature, see the attached figure. 
   >>> Set the SETTINGS
   >>> To generate 3D mesh: gmsh pipe3DMesh.geo -3 -order 2   
   >>> To generate 2D mesh: gmsh pipe3DMesh.geo -2 -order 2   
*/

//constants
PI=3.14159265359;

////////////////////////////////////////////////////////
// GENERAL SETTING /////////////////////////////////////
//***** Choose the mesh dimension
meshDim={meshDim};  //2 (2D mesh), 3 (3D mesh)
"""

    # Grid Settings Section
    grid_settings = f"""
// GRID SETTINGS ///////////////////////////////////////
//***** Geometrical parameters
// Note: r*<RB*<R
R={R};   //Pipe radius
rt={rt};
rb={rb};
ra= (rt + rb)/2.0;
RBt={RBt};
RBb={RBb};
tht={tht};  //theta top
thb={thb};  //theta bottom
lambda1t={lambda1t};   //=R_{{arc}}/R
lambda1b={lambda1b};   //=R_{{arc}}/R
lambda2={lambda2};     //=R_{{arc}}/R
dyc = {dyc};
Lz={Lz};   //length in z-dir (axial)
//***** Grid Paramaters
Nch={Nch};  // no. of nodes (=#elem+1) in azimuthal direction    # 12 16
Ncv={Ncv};  // no. of nodes (=#elem+1) in azimuthal direction    # 12 16
NB={NB};   // no. of elemtns adjacent to the wall
NM={NM};   // no. of nodes (=#elem+1) between the near wall layer and central square part # 5 7
Nc2={int((Nch+1)/2.0)}; // NM=8 for old version of mesh in gmsh
// compression ratios over the radial lines of the mesh
compressRatio_B={compressRatio_B};  //ratio of grid compression toward the wall (<1)
compressRatio_M={compressRatio_M};  //compression ratio in the middle layer
Nz={Nz};    //no of elements in z-dire (axial)
///////////////////////////////////////////////////
"""
    # Define points coordinates (top)
    geometry_creation = f"""
dxt=rt*Cos(tht);
dyt=rt*Sin(tht);
dxBt=RBt*Cos(tht);
dyBt=RBt*Sin(tht);
Dxt=R*Cos(tht);
Dyt=R*Sin(tht);
Dyxt=Hypot(dyt + lambda1t*R, dxt) - lambda1t*R;
RBC =Hypot(dyc + dyBt, dxBt) - dyc;"""
    geometry_creation += f"""
dxb=rb*Cos(thb);
dyb=rb*Sin(thb);
dxBb=RBb*Cos(thb);
dyBb=RBb*Sin(thb);
Dxb=R*Cos(thb);
Dyb=R*Sin(thb);
Dyxb=Hypot(dyb + lambda1b*R, dxb) - lambda1b*R;"""

    ipts = 0
    header = f"""
//***** define point coordinates"""
    points = [ header ]
    points.append(f"""
// inner ring right, clockwise from top""")
    ipts += 1; points.append(create_point(ipts,  0.0 , 'Dyxt', 0.0, 1.0))
    ipts += 1; points.append(create_point(ipts, 'dxt'  , 'dyt'  , 0.0, 1.0))
    ipts += 1; points.append(create_point(ipts, 'dxb'  , '-dyb' , 0.0, 1.0))
    ipts += 1; points.append(create_point(ipts, 0.0 , '-Dyxb' , 0.0, 1.0))
    points.append(f"""
// middle ring right, clockwise from top""")
    ipts += 1; points.append(create_point(ipts, 0.0, 'RBC' , 0.0, 1.0))
    ipts += 1; points.append(create_point(ipts, 'dxBt' , 'dyBt' , 0.0, 1.0))
    ipts += 1; points.append(create_point(ipts, 'dxBb' , '-dyBb', 0.0, 1.0))
    ipts += 1; points.append(create_point(ipts, 0.0, '-RBb', 0.0, 1.0))
    points.append(f"""
// outer ring right, clockwise from top""")
    ipts += 1; points.append(create_point(ipts, 0.0, ' R' , 0.0, 1.0))
    ipts += 1; points.append(create_point(ipts, 'Dxt'  , 'Dyt'  , 0.0, 1.0))
    ipts += 1; points.append(create_point(ipts, 'Dxb'  , '-Dyb' , 0.0, 1.0))
    ipts += 1; points.append(create_point(ipts, 0.0, '-R' , 0.0, 1.0))
    points.append(f"""
//auxiliary points (only help define the geometry)
// center""")
    ipts += 1; points.append(create_point(ipts,             0,             0, 0, 1.0))
    ipts += 1; points.append(create_point(ipts, pts_aux[0][0], pts_aux[0][1], 0, 1.0))
    ipts += 1; points.append(create_point(ipts,             0, '-lambda1t*R', 0, 1.0))
    ipts += 1; points.append(create_point(ipts, pts_aux[1][0], pts_aux[1][1], 0, 1.0))
    ipts += 1; points.append(create_point(ipts,             0,  'lambda1b*R', 0, 1.0))
    ipts += 1; points.append(create_point(ipts,             0,        '-dyc', 0, 1.0))
    ipts += 1; points.append(create_point(ipts, pts_aux[2][0], pts_aux[2][1], 0, 1.0))
    ipts += 1; points.append(create_point(ipts, pts_aux[3][0], pts_aux[3][1], 0, 1.0))
    if not half:
        header = f"""
// left side points"""
        points_sym = [ header ]
        ipts += 1; points_sym.append(create_point(ipts, '-dxb' , '-dyb' , 0.0, 1.0))
        ipts += 1; points_sym.append(create_point(ipts, '-dxt' , 'dyt'  , 0.0, 1.0))
        ipts += 1; points_sym.append(create_point(ipts, '-dxBb', '-dyBb', 0.0, 1.0))
        ipts += 1; points_sym.append(create_point(ipts, '-dxBt', 'dyBt' , 0.0, 1.0))
        ipts += 1; points_sym.append(create_point(ipts, '-Dxb' , '-Dyb' , 0.0, 1.0))
        ipts += 1; points_sym.append(create_point(ipts, '-Dxt' , ' Dyt' , 0.0, 1.0))
       
    # Creating circles dynamically
    icl = 0
    header = f"""
//***** define lines and curves"""
    circles = [ header ]
    circles.append(f"""
// inner ring right clockwise from top""")
    icl += 1; circles.append(create_circle(icl, 1, 15, 2))
    icl += 1; circles.append(create_circle(icl, 2, 16, 3))
    icl += 1; circles.append(create_circle(icl, 3, 17, 4))
    circles.append(f"""
// middle ring right clockwise from top""")
    icl += 1; circles.append(create_circle(icl, 5, 18, 6))
    icl += 1; circles.append(create_circle(icl, 6, 20, 7))
    icl += 1; circles.append(create_circle(icl, 7, 13, 8))
    circles.append(f"""
// outer ring right clockwise from top""")
    icl += 1; circles.append(create_circle(icl, 9, 13,10))
    icl += 1; circles.append(create_circle(icl,10, 13,11))
    icl += 1; circles.append(create_circle(icl,11, 13,12))
    
    # Creating lines dynamically
    header = f"""

// // lines in middle segment clockwise from top going outward"""
    lines = [ header ]
    icl += 1; lines.append(create_line(icl,  1,  5 ))
    icl += 1; lines.append(create_line(icl,  2,  6 ))
    icl += 1; lines.append(create_line(icl,  3,  7 ))
    icl += 1; lines.append(create_line(icl,  4,  8 ))
    lines.append(f"""
// lines in outer segment clockwise from top right""")
    icl += 1; lines.append(create_line(icl,  5,  9 ))
    icl += 1; lines.append(create_line(icl,  6, 10 ))
    icl += 1; lines.append(create_line(icl,  7, 11 ))
    icl += 1; lines.append(create_line(icl,  8, 12 ))
    lines.append(f"""
// central line upward""")
    icl += 1; lines.append(create_line(icl,  4,  1 ))
    
    if not half:
        header = f"""
// left side circles"""
        circles_sym = [ header ]
        icl += 1; circles_sym.append(create_circle(icl, 4, 17, 21))
        icl += 1; circles_sym.append(create_circle(icl,21, 14, 22))
        icl += 1; circles_sym.append(create_circle(icl,22, 15,  1))
        icl += 1; circles_sym.append(create_circle(icl, 8, 13, 23))
        icl += 1; circles_sym.append(create_circle(icl,23, 19, 24))
        icl += 1; circles_sym.append(create_circle(icl,24, 18,  5))
        icl += 1; circles_sym.append(create_circle(icl,12, 13, 25))
        icl += 1; circles_sym.append(create_circle(icl,25, 13, 26))
        icl += 1; circles_sym.append(create_circle(icl,26, 13,  9))

        header = f"""
// left side lines"""
        lines_sym = [ header ]
        icl += 1; lines_sym.append(create_line(icl, 21, 23 ))
        icl += 1; lines_sym.append(create_line(icl, 22, 24 ))
        icl += 1; lines_sym.append(create_line(icl, 23, 25 ))
        icl += 1; lines_sym.append(create_line(icl, 24, 26 ))
    
    # Creating Transfinite Lines dynamically
    header = f"""
//***** assign number of mesh on the created lines/arcs"""
    transfinite_lines = [
        header,
        create_transfinite_line([ 7, 4, 1, 3, 6, 9], 'Nc2'),
        create_transfinite_line([-18, 2, 5, 8],   'Ncv'), #, progression='compressRatio_M'),
        create_transfinite_line([10, 11, 12, 13], 'NM', progression='compressRatio_M'),
        create_transfinite_line([14, 15, 16, 17], 'NB', progression='compressRatio_B')
        ]
    
    if not half:
        transfinite_lines += [
            f"""
//*left side""",
            create_transfinite_line([ 27, 24, 21, 19, 22, 25 ], 'Nc2'),
            create_transfinite_line([-20, -23, -26], 'Ncv'), #, progression='compressRatio_M'),
            create_transfinite_line([ 28, 29 ], 'NM', progression='compressRatio_M'),
            create_transfinite_line([ 30, 31 ], 'NB', progression='compressRatio_B')
        ]
    
    # Creating Line Loops and Plane Surfaces dynamically
    header = f"""
//***** create surfaces
// Note: use a negative sign if a line is swept in the opposite direction of the original definition"""
    isf = 0
    line_loops_surfaces = [ header ]
    line_loops_surfaces.append(f"""
// central block""")
    isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [1, 2, 3, 18], isf))
    line_loops_surfaces.append(f"""
// middle ring""")
    isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [4,-11, -1, 10], isf))
    isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [5,-12, -2, 11], isf))
    isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [6,-13, -3, 12], isf))
    line_loops_surfaces.append(f"""
// outer ring""")
    isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [7,-15, -4, 14], isf))
    isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [8,-16, -5, 15], isf))
    isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [9,-17, -6, 16], isf))

    if not half:
        header = f"""

// left side"""
        line_loops_surfaces_sym = [ header ]
        isf += 1; line_loops_surfaces_sym.append(create_line_loop_surface(isf, [-18, 19, 20, 21 ], isf))
        isf += 1; line_loops_surfaces_sym.append(create_line_loop_surface(isf, [ 13, 22,-28,-19 ], isf))
        isf += 1; line_loops_surfaces_sym.append(create_line_loop_surface(isf, [ 28, 23,-29,-20 ], isf))
        isf += 1; line_loops_surfaces_sym.append(create_line_loop_surface(isf, [ 29, 24,-10,-21 ], isf))
        isf += 1; line_loops_surfaces_sym.append(create_line_loop_surface(isf, [ 17, 25,-30,-22 ], isf))
        isf += 1; line_loops_surfaces_sym.append(create_line_loop_surface(isf, [ 30, 26,-31,-23 ], isf))
        isf += 1; line_loops_surfaces_sym.append(create_line_loop_surface(isf, [ 31, 27,-14,-24 ], isf))

    # Combine all circles, lines, transfinite lines, and line loops into the script
    geometry_creation += "\n".join(points) + "\n"
    if not half:
        geometry_creation += "\n".join(points_sym) + "\n"
    geometry_creation += "\n".join(circles) + "\n"
    geometry_creation += "\n".join(lines) + "\n"
    if not half:
        geometry_creation += "\n".join(circles_sym) + "\n"
        geometry_creation += "\n".join(lines_sym) + "\n"
    geometry_creation += "\n".join(transfinite_lines) + "\n"
    geometry_creation += "\n".join(line_loops_surfaces) + "\n"
    if not half:
        geometry_creation += "\n".join(line_loops_surfaces_sym) + "\n"

    if half:
       case2d = f"""
If (meshDim==2)
   Physical Line("wall")={{7, 8, 9}};
   Physical Line("sym")={{14, 10, 18, 13, 17}};
   Physical Surface(1)={{1:7}};
EndIf"""
    else:
        case2d = f"""
If (meshDim==2)
   Physical Line("wall")={{7, 8, 9, 25, 26, 27}};
   Physical Surface(1)={{1:14}};
EndIf"""
    
    surface = f"""
Recombine Surface "*";
Transfinite Surface "*";"""
    
    case3d = f"""
If (meshDim==3)

   //make a 3d mesh by extrusion in z-dir
   mesh3D[]=Extrude {{0,0,Lz}} 
   {{
       Surface{{1:9}};
       Layers{{Nz}}; 
       Recombine; 
   }};

   //Physical Surfaces & Volume (Note: gmsh only generates mesh for the physical entities)
   // BC tag of the surfaces are assigned in accordance with what is added in usrdat2() routine in case.usr. This is in accordance with the requirements by gmsh2nek. see the following link:
   //https://github.com/yhaomin2007/Nek5000/tree/master/gmsh2nek_sourcecode/gmsh2nek/
   // 1: inlet
   // 2: outlet
   // 3: wall
   Physical Surface("inlet") = {{6, 7, 8, 9, 5, 2, 3, 4, 1}};
   Physical Surface("outlet") = {{152, 174, 196, 218, 64, 86, 108, 130, 42}};
   Physical Surface("wall") = {{139, 213, 191, 165}};
   Physical Volume("flowDomain") = {{6, 2, 1, 4, 8, 3, 7, 5, 9}};

   Recombine Volume "*";

EndIf
    """

    savemesh = f"""
Coherence;

/////////////////////////////////////////////////////////////////////
// Mesh saving section
////////////////////////////////////////////////////////////////////
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 0;
Mesh.Binary = 0;"""

    # Write the script to a file
    with open(filename, "w") as f:
        f.write(general_settings)
        f.write(grid_settings)
        f.write(geometry_creation)
        f.write(case2d)
        f.write(surface)
        f.write(case3d)
        f.write(savemesh)

    print(f"Mesh script saved to {filename}")
