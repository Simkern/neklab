
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
rt=0.5;
rb=0.7;
ra= (rt + rb)/2.0;
RBt=0.92;
RBb=0.94;
tht=0.5235987755982988;  //theta top
thb=1.0471975511965976;  //theta bottom
lambda1t=0.75;   //=R_{arc}/R
lambda1b=0.6;   //=R_{arc}/R
lambda2=0.5;     //=R_{arc}/R
dyc = 0.04;
Lz=1;   //length in z-dir (axial)
//***** Grid Paramaters
Nch=11;  // no. of nodes (=#elem+1) in azimuthal direction    # 12 16
Ncv=11;  // no. of nodes (=#elem+1) in azimuthal direction    # 12 16
NB=1;   // no. of elemtns adjacent to the wall
NM=5;   // no. of nodes (=#elem+1) between the near wall layer and central square part # 5 7
Nc2=6;
     // NM=8 for old version of mesh in gmsh
// compression ratios over the radial lines of the mesh
compressRatio_B=0.6;  //ratio of grid compression toward the wall (<1)
compressRatio_M=0.95;  //compression ratio in the middle layer
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
//***** define points coordinates
//auxiliary points (only help define the geometry)
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.6161914625146543, -0.08040977798980793, 0, 1.0};
Point(3) = {0, -lambda1t*R, 0, 1.0};
Point(4) = {-0.6161914625146543, -0.08040977798980793, 0, 1.0};
Point(5) = {0, lambda1b*R, 0, 1.0};
Point(6) = {0, -dyc, 0, 1.0};
Point(7) = {0.017214092452766394, -0.010184269096062515, 0, 1.0};
Point(8) = {-0.017214092452766394, -0.010184269096062515, 0, 1.0};

//blocks vertices
Point(9) = {dxt, dyt, 0.0, 1.0};
Point(10) = {dxb, -dyb, 0.0, 1.0};
Point(11) = {0.0, -Dyxb, 0.0, 1.0};
Point(12) = {0.0, Dyxt, 0.0, 1.0};
Point(13) = {dxBt, dyBt, 0.0, 1.0};
Point(14) = {dxBb, -dyBb, 0.0, 1.0};
Point(15) = {0.0, -RBb, 0.0, 1.0};
Point(16) = {0.0, RBC, 0.0, 1.0};
Point(17) = {Dxt, Dyt, 0.0, 1.0};
Point(18) = {Dxb, -Dyb, 0.0, 1.0};
Point(19) = {0.0, -R, 0.0, 1.0};
Point(20) = {0.0,  R, 0.0, 1.0};

//***** define lines and curves
Circle(1)={12, 3, 9};
Circle(2)={9, 4, 10};
Circle(3)={10, 5, 11};
Circle(4)={16, 6, 13};
Circle(5)={13, 8, 14};
Circle(6)={14, 1, 15};
Circle(7)={20, 1, 17};
Circle(8)={17, 1, 18};
Circle(9)={18, 1, 19};
Line(10)={ 9, 13 };
Line(11)={ 10, 14 };
Line(12)={ 11, 15 };
Line(13)={ 12, 16 };
Line(14)={ 13, 17 };
Line(15)={ 14, 18 };
Line(16)={ 15, 19 };
Line(17)={ 16, 20 };
Line(18)={ 11, 12 };

//***** assign number of mesh on the created lines/arcs
Transfinite Line { 7, 4, 1, 3, 6, 9 } = Nc2;
Transfinite Line { -18, 2, 5, 8 } = Ncv Using Progression compressRatio_M;
Transfinite Line { 10, 11, 12, 13 } = NM Using Progression compressRatio_M;
Transfinite Line { 14, 15, 16, 17 } = NB Using Progression compressRatio_B;

//***** create surfaces
// Note: use a negative sign if a line is swept in the opposite direction of the original definition
Line Loop(1)={ 1, 2, 3, 18 };   Plane Surface(1)={ 1 };
Line Loop(2)={ 4, -10, -1, 13 };   Plane Surface(2)={ 2 };
Line Loop(3)={ 7, -14, -4, 17 };   Plane Surface(3)={ 3 };
Line Loop(4)={ 2, 11, -5, -10 };   Plane Surface(4)={ 4 };
Line Loop(5)={ 5, 15, -8, -14 };   Plane Surface(5)={ 5 };
Line Loop(6)={ 3, 12, -6, -11 };   Plane Surface(6)={ 6 };
Line Loop(7)={ 6, 16, -9, -15 };   Plane Surface(7)={ 7 };

If (meshDim==2)
   Physical Line("wall")={7, 8, 9};
   Physical Line("sym")={16, 12, 18, 13, 17};
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