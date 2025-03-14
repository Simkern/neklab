
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
meshDim=2;  //2 (2D mesh), 3 (3D mesh)

// GRID SETTINGS ///////////////////////////////////////
//***** Geometrical parameters
// Note: r*<RB*<R
R=1.0;   //Pipe radius
rt=0.7;
rb=0.7;
ra= (rt + rb)/2.0;
RBt=0.92;
RBb=0.92;
tht=0.7853981633974483;  //theta top
thb=0.7853981633974483;  //theta bottom
lambda1t=0.8;   //=R_{arc}/R
lambda1b=0.8;   //=R_{arc}/R
lambda2=0.8;     //=R_{arc}/R
dyc = 0.0;
Lz=1;   //length in z-dir (axial)
//***** Grid Paramaters
Nch=7;  // no. of nodes (=#elem+1) in azimuthal direction    # 12 16
Ncv=7;  // no. of nodes (=#elem+1) in azimuthal direction    # 12 16
NB=1;   // no. of elemtns adjacent to the wall
NM=3;   // no. of nodes (=#elem+1) between the near wall layer and central square part # 5 7
Nc2=4; // NM=8 for old version of mesh in gmsh
// compression ratios over the radial lines of the mesh
compressRatio_B=0.85;  //ratio of grid compression toward the wall (<1)
compressRatio_M=0.87;  //compression ratio in the middle layer
Nz=180;    //no of elements in z-dire (axial)
///////////////////////////////////////////////////

dxt=rt*Cos(tht);
dyt=rt*Sin(tht);
dxBt=RBt*Cos(tht);
dyBt=RBt*Sin(tht);
Dxt=R*Cos(tht);
Dyt=R*Sin(tht);
Dyxt=Hypot(dyt + lambda1t*R, dxt) - lambda1t*R;
RBC =Hypot(dyc + dyBt, dxBt) - dyc;
dxb=rb*Cos(thb);
dyb=rb*Sin(thb);
dxBb=RBb*Cos(thb);
dyBb=RBb*Sin(thb);
Dxb=R*Cos(thb);
Dyb=R*Sin(thb);
Dyxb=Hypot(dyb + lambda1b*R, dxb) - lambda1b*R;
//***** define point coordinates

// inner ring right, clockwise from top
Point(1) = {0.0, Dyxt, 0.0, 1.0};
Point(2) = {dxt, dyt, 0.0, 1.0};
Point(3) = {dxb, -dyb, 0.0, 1.0};
Point(4) = {0.0, -Dyxb, 0.0, 1.0};

// middle ring right, clockwise from top
Point(5) = {0.0, RBC, 0.0, 1.0};
Point(6) = {dxBt, dyBt, 0.0, 1.0};
Point(7) = {dxBb, -dyBb, 0.0, 1.0};
Point(8) = {0.0, -RBb, 0.0, 1.0};

// outer ring right, clockwise from top
Point(9) = {0.0,  R, 0.0, 1.0};
Point(10) = {Dxt, Dyt, 0.0, 1.0};
Point(11) = {Dxb, -Dyb, 0.0, 1.0};
Point(12) = {0.0, -R, 0.0, 1.0};

//auxiliary points (only help define the geometry)
// center
Point(13) = {0, 0, 0, 1.0};
Point(14) = {0.9210054790200463, 0.0, 0, 1.0};
Point(15) = {0, -lambda1t*R, 0, 1.0};
Point(16) = {-0.9210054790200463, 0.0, 0, 1.0};
Point(17) = {0, lambda1b*R, 0, 1.0};
Point(18) = {0, -dyc, 0, 1.0};
Point(19) = {1.1102230246251565e-16, 0.0, 0, 1.0};
Point(20) = {-1.1102230246251565e-16, 0.0, 0, 1.0};

//***** define lines and curves

// inner ring right clockwise from top
Circle(1)={1, 15, 2};
Circle(2)={2, 16, 3};
Circle(3)={3, 17, 4};

// middle ring right clockwise from top
Circle(4)={5, 18, 6};
Circle(5)={6, 20, 7};
Circle(6)={7, 13, 8};

// outer ring right clockwise from top
Circle(7)={9, 13, 10};
Circle(8)={10, 13, 11};
Circle(9)={11, 13, 12};


// // lines in middle segment clockwise from top going outward
Line(10)={ 1, 5 };
Line(11)={ 2, 6 };
Line(12)={ 3, 7 };
Line(13)={ 4, 8 };

// lines in outer segment clockwise from top right
Line(14)={ 5, 9 };
Line(15)={ 6, 10 };
Line(16)={ 7, 11 };
Line(17)={ 8, 12 };

// central line upward
Line(18)={ 4, 1 };

//***** assign number of mesh on the created lines/arcs
Transfinite Line { 7, 4, 1, 3, 6, 9 } = Nc2;
Transfinite Line { -18, 2, 5, 8 } = Ncv Using Progression compressRatio_M;
Transfinite Line { 10, 11, 12, 13 } = NM Using Progression compressRatio_M;
Transfinite Line { 14, 15, 16, 17 } = NB Using Progression compressRatio_B;

//***** create surfaces
// Note: use a negative sign if a line is swept in the opposite direction of the original definition

// central block
Line Loop(1)={ 1, 2, 3, 18 };   Plane Surface(1)={ 1 };

// middle ring
Line Loop(2)={ 4, -11, -1, 10 };   Plane Surface(2)={ 2 };
Line Loop(3)={ 5, -12, -2, 11 };   Plane Surface(3)={ 3 };
Line Loop(4)={ 6, -13, -3, 12 };   Plane Surface(4)={ 4 };

// outer ring
Line Loop(5)={ 7, -15, -4, 14 };   Plane Surface(5)={ 5 };
Line Loop(6)={ 8, -16, -5, 15 };   Plane Surface(6)={ 6 };
Line Loop(7)={ 9, -17, -6, 16 };   Plane Surface(7)={ 7 };

If (meshDim==2)
   Physical Line("wall")={7, 8, 9};
   Physical Line("sym")={14, 10, 18, 13, 17};
   Physical Surface(1)={1:7};
EndIf
Recombine Surface "*";
Transfinite Surface "*";
If (meshDim==3)

   //make a 3d mesh by extrusion in z-dir
   mesh3D[]=Extrude {0,0,Lz} 
   {
       Surface{1:9};
       Layers{Nz}; 
       Recombine; 
   };

   //Physical Surfaces & Volume (Note: gmsh only generates mesh for the physical entities)
   // BC tag of the surfaces are assigned in accordance with what is added in usrdat2() routine in case.usr. This is in accordance with the requirements by gmsh2nek. see the following link:
   //https://github.com/yhaomin2007/Nek5000/tree/master/gmsh2nek_sourcecode/gmsh2nek/
   // 1: inlet
   // 2: outlet
   // 3: wall
   Physical Surface("inlet") = {6, 7, 8, 9, 5, 2, 3, 4, 1};
   Physical Surface("outlet") = {152, 174, 196, 218, 64, 86, 108, 130, 42};
   Physical Surface("wall") = {139, 213, 191, 165};
   Physical Volume("flowDomain") = {6, 2, 1, 4, 8, 3, 7, 5, 9};

   Recombine Volume "*";

EndIf
    
Coherence;

/////////////////////////////////////////////////////////////////////
// Mesh saving section
////////////////////////////////////////////////////////////////////
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 0;
Mesh.Binary = 0;