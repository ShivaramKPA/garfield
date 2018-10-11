# garfield

To compile GMSH, Elmer, and Garfield:

- RTPC.geo which is the file used to define the geometry for the finite 
element method (FEM). 

- RTPC.msh which is the mesh file. The mesh is the list of all finite 
elements in our geometry. It has been created with:

> gmsh RTPC.geo -3 -order 2

- RTPC/mesh.boundary, RTPC/mesh.elements, RTPC/mesh.header, and 
RTPC/mesh.nodes. These files will be read by ElmerSolver to create the 
electric file and later by Garfield. They were created with:

> ElmerGrid 14 2 <mesh_file>.msh -autoclean

- RTPC.sif, which is the file defining the physics to apply to the 
mesh. The boundary conditions (values of potentials) are implemented 
there.

- RTPC/RTPC.result and RTPC/RTPC.ep. Only the .result file will be 
relevant for simulations with Garfield, as the .ep file is meant to be 
read into a visualization program. They were created by:

> ElmerSolver RTPC.sif

which performs the calculation described by the given solver input file 
(sif).

- RTPC.C and makefile which are the files used for the Garfield 
simulation. To run it do:

> make 

(it will create a executable called RTPC_garf)

> ./garf_RTPC
