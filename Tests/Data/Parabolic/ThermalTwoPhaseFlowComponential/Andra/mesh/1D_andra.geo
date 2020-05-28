// Gmsh project created on Wed Aug 07 14:36:13 2019
lc=0.05;
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {7, 0, 0, lc};
//+
Point(3) = {10, 0, 0, lc};


//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
Transfinite Line{1} = 141;
Transfinite Line{2} = 300 Using Progression 0.98;
//+
Physical Line(1) = {1,2};
