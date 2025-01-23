import sys, os, re
import copy
import numpy as np
from itertools import product
import matplotlib.pyplot as plt

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

def generate_gmsh_script(R, rt, rb, RBt, RBb, tht, thb, lambda1, lambda2, dyc, Lz, Nch, Ncv, NB, NM, compressRatio_B, compressRatio_M, Nz, half=False, meshDim=2, filename="pipe_mesh.geo"):

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

    ## INNER RING
    # midpoint
    xm = (dxt + dxb)/2.0
    ym = (dyt - dyb)/2.0
    # bisector
    dx = dxb - dxt
    dy = - dyb - dyt
    #m  = -dx/dy
    norm = np.sqrt(dx**2 + dy**2)
    h   = norm/2.0
    ux  = dy/norm
    uy  = -dx/norm
    # distance along bisector
    L   = np.sqrt((lambda2*R + ra)**2 - h**2)

    paux0 = [ -xm - L*ux, ym + L*uy ]
    paux1 = [  xm + L*ux, ym + L*uy ]
    
    ## MIDDLE RING
    # midpoint
    xmB = (dxBt + dxBb)/2.0
    ymB = (dyBt - dyBb)/2.0
    # bisector
    dxB = dxBb - dxBt
    dyB = - dyBb - dyBt
    #m  = -dx/dy
    norm = np.sqrt(dxB**2 + dyB**2)
    hB   = norm/2.0
    uxB  = dyB/norm
    uyB  = -dxB/norm
    # distance along bisector
    LB   = np.sqrt(RBb**2 - hB**2)

    paux2 = [ -xmB - LB*uxB, ymB + LB*uyB ]
    paux3 = [  xmB + LB*uxB, ymB + LB*uyB ]

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
// Note: r<RB<R
R={R};   //Pipe radius
rt={rt};
rb={rb};
ra= (rt + rb)/2.0;
RBt={RBt};    // 0.95                                                 //0.97;   
RBb={RBb};    // 0.95                                                 //0.97;   
tht={tht};  //theta top
thb={thb};  //theta bottom
lambda1={lambda1};   //=R_{{arc}}/R   0.3
lambda2={lambda2};   //=R_{{arc}}/R   0.3
dyc = {dyc};
Lz={Lz};   //length in z-dir (axial)
//***** Grid Paramaters
Nch={Nch};  // no. of nodes (=#elem+1) in azimuthal direction    # 12 16
Ncv={Ncv};  // no. of nodes (=#elem+1) in azimuthal direction    # 12 16
NB={NB};   // no. of elemtns adjacent to the wall
NM={NM};   // no. of nodes (=#elem+1) between the near wall layer and central square part # 5 7
Nc2={int((Nch+1)/2.0)};
     // NM=8 for old version of mesh in gmsh
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
Dyt=R*Sin(tht);"""
    if (half):
        geometry_creation += f"""
Dyxt=Hypot(dyt + lambda1*R, dxt) - lambda1*R;
RBC =Hypot(dyc + dyBt, dxBt) - dyc;"""
    geometry_creation += f"""
dxb=rb*Cos(thb);
dyb=rb*Sin(thb);
dxBb=RBb*Cos(thb);
dyBb=RBb*Sin(thb);
Dxb=R*Cos(thb);
Dyb=R*Sin(thb);"""
    if (half):
        geometry_creation += f"""
Dyxb=Hypot(dyb + lambda1*R, dxb) - lambda1*R;"""

    ipts = 0
    header = f"""
//***** define points coordinates
//auxiliary points (only help define the geometry)"""
    aux_points = [ header ]
    ipts += 1; aux_points.append(create_point(ipts,        0,          0, 0, 1.0))
         #create_point(2, 'lambda*R', 0, 0, 1.0),
    ipts += 1; aux_points.append(create_point(ipts, paux0[0],   paux0[1], 0, 1.0))
    ipts += 1; aux_points.append(create_point(ipts,        0,'-lambda1*R', 0, 1.0))
         #create_point(4, '-lambda*R', 0, 0, 1.0),
    ipts += 1; aux_points.append(create_point(ipts, paux1[0],   paux1[1], 0, 1.0))
    ipts += 1; aux_points.append(create_point(ipts,        0, 'lambda1*R', 0, 1.0))
    ipts += 1; aux_points.append(create_point(ipts,        0,     '-dyc', 0, 1.0))
    ipts += 1; aux_points.append(create_point(ipts, paux2[0],   paux2[1], 0, 1.0))
    ipts += 1; aux_points.append(create_point(ipts, paux3[0],   paux3[1], 0, 1.0))

    naux = ipts
    header = f"""
//blocks vertices"""
    block_points = [ header ]
    ipts += 1; block_points.append(create_point(ipts, 'dxt'  , 'dyt'  , 0.0, 1.0))
    ipts += 1; block_points.append(create_point(ipts, 'dxb'  , '-dyb' , 0.0, 1.0))
    if (half):
        ipts += 1; block_points.append(create_point(ipts, 0.0 , '-Dyxb' , 0.0, 1.0))
        ipts += 1; block_points.append(create_point(ipts, 0.0 , 'Dyxt'  , 0.0, 1.0))
    else:
        ipts += 1; block_points.append(create_point(ipts, '-dxb' , '-dyb' , 0.0, 1.0))
        ipts += 1; block_points.append(create_point(ipts, '-dxt' , 'dyt'  , 0.0, 1.0))
    ipts += 1; block_points.append(create_point(ipts, 'dxBt' , 'dyBt' , 0.0, 1.0))
    ipts += 1; block_points.append(create_point(ipts, 'dxBb' , '-dyBb', 0.0, 1.0))
    if (half):
        ipts += 1; block_points.append(create_point(ipts, 0.0, '-RBb', 0.0, 1.0))
        ipts += 1; block_points.append(create_point(ipts, 0.0, 'RBC' , 0.0, 1.0))
    else:
        ipts += 1; block_points.append(create_point(ipts, '-dxBb', '-dyBb', 0.0, 1.0))
        ipts += 1; block_points.append(create_point(ipts, '-dxBt', 'dyBt' , 0.0, 1.0))
    ipts += 1; block_points.append(create_point(ipts, 'Dxt'  , 'Dyt'  , 0.0, 1.0))
    ipts += 1; block_points.append(create_point(ipts, 'Dxb'  , '-Dyb' , 0.0, 1.0))
    if (half):
        ipts += 1; block_points.append(create_point(ipts, 0.0, '-R' , 0.0, 1.0))
        ipts += 1; block_points.append(create_point(ipts, 0.0, ' R' , 0.0, 1.0))
    else:
        ipts += 1; block_points.append(create_point(ipts, '-Dxb' , '-Dyb' , 0.0, 1.0))
        ipts += 1; block_points.append(create_point(ipts, '-Dxt' , ' Dyt' , 0.0, 1.0))

    # Creating circles dynamically
    icl = 0
    header = f"""
//***** define lines and curves"""
    circles = [ header ]
    icl += 1;     circles.append(create_circle(icl, naux+4 , 3, naux+1))
    icl += 1;     circles.append(create_circle(icl, naux+1 , 4, naux+2))
    icl += 1;     circles.append(create_circle(icl, naux+2 , 5, naux+3))
    if not half:
        icl += 1; circles.append(create_circle(icl, naux+3 , 2, naux+4))
    icl += 1;     circles.append(create_circle(icl, naux+8 , 6, naux+5))
    icl += 1;     circles.append(create_circle(icl, naux+5 , 8, naux+6))
    icl += 1;     circles.append(create_circle(icl, naux+6 , 1, naux+7))
    if not half:
        icl += 1; circles.append(create_circle(icl, naux+7 , 7, naux+8))
    icl += 1;     circles.append(create_circle(icl, naux+12, 1, naux+9))
    icl += 1;     circles.append(create_circle(icl, naux+9 , 1, naux+10))
    icl += 1;     circles.append(create_circle(icl, naux+10, 1, naux+11))
    if not half:
        icl += 1; circles.append(create_circle(icl, naux+11, 1, naux+12))
    
    # Creating lines dynamically
    lines = [ ]
    icl += 1;     lines.append(create_line(icl, naux+1, naux+5))
    icl += 1;     lines.append(create_line(icl, naux+2, naux+6))
    icl += 1;     lines.append(create_line(icl, naux+3, naux+7))
    icl += 1;     lines.append(create_line(icl, naux+4, naux+8))
    icl += 1;     lines.append(create_line(icl, naux+5, naux+9))
    icl += 1;     lines.append(create_line(icl, naux+6, naux+10))
    icl += 1;     lines.append(create_line(icl, naux+7, naux+11))
    icl += 1;     lines.append(create_line(icl, naux+8, naux+12))
    if (half):
        icl += 1; lines.append(create_line(icl, naux+3, naux+4))
    
    # Creating Transfinite Lines dynamically
    header = f"""
//***** assign number of mesh on the created lines/arcs"""
    if (half):
        transfinite_lines = [
            header,
            #create_transfinite_line([1, 4, 7], 'Nc2'),
            #create_transfinite_line([2, 5, 8], 'Nc2'),
            #create_transfinite_line([3, 6, 9], 'Nc2'),
            #create_transfinite_line([18, 2, 5, 8], 'Nc'),
            create_transfinite_line([ 7, 4, 1, 3, 6, 9], 'Nc2'),
            create_transfinite_line([-18, 2, 5, 8],   'Ncv', progression='compressRatio_M'),
            create_transfinite_line([10, 11, 12, 13], 'NM', progression='compressRatio_M'),
            create_transfinite_line([14, 15, 16, 17], 'NB', progression='compressRatio_B')
        ]
    else:
        transfinite_lines = [
            header,
            #create_transfinite_line([1, 2, 3, 4], 'Nc', bump=1.0),
            #create_transfinite_line([5, 6, 7, 8], 'Nc'),
            #create_transfinite_line([9, 10, 11, 12], 'Nc'),
            create_transfinite_line([ 9, 5, 1, 3, 7, 11], 'Nch'),
            create_transfinite_line([-12, -8, -4, 2, 6, 10], 'Ncv', progression='compressRatio_M'),
            create_transfinite_line([13, 14, 15, 16], 'NM', progression='compressRatio_M'),
            create_transfinite_line([17, 18, 19, 20], 'NB', progression='compressRatio_B')
        ]
    
    # Creating Line Loops and Plane Surfaces dynamically
    header = f"""
//***** create surfaces
// Note: use a negative sign if a line is swept in the opposite direction of the original definition"""
    isf = 0
    if (half):
        line_loops_surfaces = [ header ]
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [1, 2, 3, 18], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [4, -10, -1, 13], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [7, -14, -4, 17], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [2, 11, -5, -10], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [5, 15, -8, -14], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [3, 12, -6, -11], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [6, 16, -9, -15], isf))
    else:
        line_loops_surfaces = [ header ]
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [1, 2, 3, 4], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [5, -13, -1, 16], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [13, 6, -14, -2], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [-3, 14, 7, -15], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [-16, -4, 15, 8], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [9, -17, -5, 20], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [17, 10, -18, -6], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [-7, 18, 11, -19], isf))
        isf += 1; line_loops_surfaces.append(create_line_loop_surface(isf, [-8, 19, 12, -20], isf))

    # Combine all circles, lines, transfinite lines, and line loops into the script
    geometry_creation += "\n".join(aux_points) + "\n"
    geometry_creation += "\n".join(block_points) + "\n"
    geometry_creation += "\n".join(circles) + "\n"
    geometry_creation += "\n".join(lines) + "\n"
    geometry_creation += "\n".join(transfinite_lines) + "\n"
    geometry_creation += "\n".join(line_loops_surfaces) + "\n"

    if (half):
       case2d = f"""
If (meshDim==2)
   Physical Line("wall")={{7, 8, 9}};
   Physical Line("sym")={{16, 12, 18, 13, 17}};
   Physical Surface(1)={{1:7}};
EndIf"""
    else:
        case2d = f"""
If (meshDim==2)
   Physical Line("wall")={{9, 10, 11, 12}};
   Physical Surface(1)={{1:9}};
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

if __name__ == "__main__":
   # Given parameters
   half = False
   R = 1.0
   r = 0.7
   RB = 0.92
   th = np.pi / 4.0
   #thh = th / 2.0
   lambda_val = 0.8
   Lz = 1
   Nc = 7
   NB = 1
   NM = 3
   compressRatio_B = 0.85
   compressRatio_M = 0.87
   Nz = 180

   # Calculate coordinates based on the formulas
   dx = r * np.cos(th)
   dy = r * np.sin(th)
   dxB = RB * np.cos(th)
   dyB = RB * np.sin(th)
   Dx = R * np.cos(th)
   Dy = R * np.sin(th)
   Dyx = np.sqrt((dx + lambda_val*R)**2 + dy**2) - lambda_val*R

   points_aux = np.array([
            [0, 0],
            [lambda_val * R, 0],
            [0, -lambda_val * R],
            [-lambda_val * R, 0],
            [0, lambda_val * R]
   ])

   # Block vertices
   if (half):
    points_block = np.array([
        [dx, dy],
        [dx, -dy],
        [0.0, -Dyx],
        [0.0, Dyx],
        [dxB, dyB],
        [dxB, -dyB],
        [0.0, -RB],
        [0.0, RB],
        [Dx, Dy],
        [Dx, -Dy],
        [0.0, -R],
        [0.0, R]
    ])
   else:
    points_block = np.array([
      [dx, dy],
      [dx, -dy],
      [-dx, -dy],
      [-dx, dy],
      [dxB, dyB],
      [dxB, -dyB],
      [-dxB, -dyB],
      [-dxB, dyB],
      [Dx, Dy],
      [Dx, -Dy],
      [-Dx, -Dy],
      [-Dx, Dy]
   ])

   points = np.concatenate([points_aux, points_block], axis=0)

   # Circle connections (correspond to the given circle definitions)
   if (half):
    circles = [
        [9, 3, 6],
        [6, 4, 7],
        [7, 5, 8],
        [13, 1, 10],
        [10, 1, 11],
        [11, 1, 12],
        [17, 1, 14],
        [14, 1, 15],
        [15, 1, 16]
    ]
   else:
      circles = [
        [9, 3, 6],
        [6, 4, 7],
        [7, 5, 8],
        [9, 2, 8],
        [13, 1, 10],
        [10, 1, 11],
        [11, 1, 12],
        [13, 1, 12],
        [17, 1, 14],
        [14, 1, 15],
        [15, 1, 16],
        [17, 1, 16]
    ]

   # Line connections (correspond to the given line definitions)
   lines = [
      [6, 10],
      [7, 11],
      [8, 12],
      [9, 13],
      [10, 14],
      [11, 15],
      [12, 16],
      [13, 17],
   ]
   if (half):
      lines = np.concatenate([ lines, [[8, 9]] ], axis=0)

   
   # Plot the points
   fig, ax = plt.subplots(figsize=(10,8))

   # Plot auxiliary points (group 1)
   ax.scatter(points_aux[:, 0], points_aux[:, 1], color='blue', label='Auxiliary Points')

   # Plot block vertices (group 2)
   ax.scatter(points_block[:, 0], points_block[:, 1], color='red', label='Block Vertices')

   # Annotate auxiliary points
   for i, point in enumerate(points_aux):
      ax.text(point[0], point[1], f'{i+1}', color='blue', fontsize=12, ha='right', va='bottom')

   # Annotate block vertices
   for i, point in enumerate(points_block):
      ax.text(point[0], point[1], f'{points_aux.shape[0] + i+1}', color='red', fontsize=12, ha='right', va='bottom')

   # Plot circles (connecting points according to the circle definitions)
   for cidx, circle in enumerate(circles):
      p1 = points[circle[0] - 1]  # Subtract 1 for 0-indexing
      c  = points[circle[1] - 1]    #center
      p2 = points[circle[2] - 1]

      cc  = c[0] + 1j*c[1]
      p1c = p1[0] + 1j*p1[1]
      p2c = p2[0] + 1j*p2[1]
      r1  = p1c - cc
      r2  = p2c - cc
      da  = np.angle(r2) - np.angle(r1)
      if (abs(da) > np.pi):
         da = 2*np.pi - abs(da)
      alp = np.linspace(0,da,101, endpoint=True)
      circ = cc + r1*np.exp(1j*alp)
      
      # Draw the circle segments
      #ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='green', linestyle='-', linewidth=1)
      ax.plot(np.real(circ), np.imag(circ), color='green', linestyle='-', linewidth=1)
         
      # Add annotation at the midpoint
      ax.text(np.real(circ[50]), np.imag(circ[50]), f'{cidx+1}', color='green', fontsize=12, ha='left', va='bottom')

   # Plot lines (connecting points according to the line definitions)
   for lidx, line in enumerate(lines):
      p1 = points[line[0] - 1]  # Subtract 1 for 0-indexing
      p2 = points[line[1] - 1]
      
      # Draw the line connecting the two points
      ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='purple', linestyle='-', linewidth=1)

      # Calculate the midpoint for annotation
      midpoint = [(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2]
      
      # Add annotation at the midpoint
      ax.text(midpoint[0], midpoint[1], f'{cidx+1+lidx+1}', color='purple', fontsize=12, ha='left', va='bottom')

   # Example usage:
   #generate_gmsh_script(R, r, RB, 'PI/4.', lambda_val, Lz, Nc, NB, NM, compressRatio_B, compressRatio_M, half=True)

   # Labels and title
   ax.set_xlabel('X')
   ax.set_ylabel('Y')
   ax.set_title('Plot of Points')
   ax.legend()

   # Display the plot
   plt.axis('equal')
   plt.show()