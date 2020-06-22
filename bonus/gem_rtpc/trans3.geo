p1 = newp;
Point(p1) = {0.0,0,-20,0.5};
//--GEM3
p2 = newp; Point(p2) = {7.61,0,-20,0.5};
p3 = newp; Point(p3) = {0,7.61,-20,0.5};
p4 = newp; Point(p4) = {-7.61,0,-20,0.5};
p5 = newp; Point(p5) = {0,-7.61,-20,0.5};

c1 = newc; Circle(c1) = {p2,p1,p3};
c2 = newc; Circle(c2) = {p3,p1,p4};
c3 = newc; Circle(c3) = {p4,p1,p5};
c4 = newc; Circle(c4) = {p5,p1,p2};

l1 = newl; Line Loop(l1) = {c1,c2,c3,c4};
//--
//--Readout Pad Board (RPB)
p6 = newp; Point(p6) = {8.0,0,-20,0.5};
p7 = newp; Point(p7) = {0,8.0,-20,0.5};
p8 = newp; Point(p8) = {-8.0,0,-20,0.5};
p9 = newp; Point(p9) = {0,-8.0,-20,0.5};

c5 = newc; Circle(c5) = {p6,p1,p7};
c6 = newc; Circle(c6) = {p7,p1,p8};
c7 = newc; Circle(c7) = {p8,p1,p9};
c8 = newc; Circle(c8) = {p9,p1,p6};

l2 = newl; Line Loop(l2) = {c5,c6,c7,c8};
//--

//--Create the surface between the two circles
s1 = news; Plane Surface(s1) = {l1,l2};

//--Extrude this surface to create the drift volume
vol_index[] = Extrude {0,0,40}{
  Surface{s1};
};
//--

//--Create the physical entities
//--The two surfaces for the GEM3 and RPB foils
pphys1 = newreg; Physical Surface("GEM3",pphys1) = {32, 24, 28, 36};
pphys2 = newreg; Physical Surface("RPB",pphys2) = {48, 44, 52, 40};

//--The drift volume
pphys3 = newreg; Physical Volume("GEM3-RPB Drift Area",pphys3) = {vol_index[1]};
Printf("New stuff %g and %g and %g",pphys1,pphys2,pphys3);
//--
