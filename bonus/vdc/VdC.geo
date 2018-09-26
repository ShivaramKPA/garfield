lc = 0.5;

p1 = newp;
Point(p1) = {0.0,0,0,lc};
//-- Create wall
p2 = newp; Point(p2) = {1,-3,0,lc};
p3 = newp; Point(p3) = {1,3,0,lc};
p4 = newp; Point(p4) = {8,3,0,lc};
p5 = newp; Point(p5) = {8,-3,0,lc};

l1 = newl; Line(l1) = {p2,p3};
l2 = newl; Line(l2) = {p3,p4};
l3 = newl; Line(l3) = {p4,p5};
l4 = newl; Line(l4) = {p5,p2};

ll1 = newl; Line Loop(ll1) = {l1,l2,l3,l4};

//--Create the rectangular surface
s1 = news; Plane Surface(s1) = {ll1};

//--Extrude this surface to create the drift volume
vol_index[] = Extrude {0,0,8}{
  Surface{s1};
};
//--

//--Create the physical entities
//--The two surfaces for the cathode and anode foils
//--pphys1 = newreg; Physical Surface("Cathode",pphys1) = {32, 24, 28, 36};
//-- pphys2 = newreg; Physical Surface("Anode",pphys2) = {48, 44, 52, 40};

//--The drift volume
//--pphys3 = newreg; Physical Volume("Drift Area",pphys3) = {vol_index[1]};
//--Printf("New stuff %g and %g and %g",pphys1,pphys2,pphys3);
//--
