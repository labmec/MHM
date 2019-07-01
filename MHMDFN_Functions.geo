// Gmsh project created on Sun Oct 28 15:40:30 2018
SetFactory("OpenCASCADE");
//+
size = 0.02;

Point(1) = {0, 0, 0, size};
Point(2) = {1, 0, 0, size};
Point(3) = {1, 1, 0, size};
Point(4) = {0, 1, 0, size};
Point(5) = {0.5, 0, 0, size/4};
Point(6) = {0.6, 0.6, 0, size/4};

Line(1) = {1,5};
Line(2) = {5,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,1};
Line(6) = {5,6};
//+
Physical Point("BC_POINT") = {1};
Physical Point("BC_FRACPOINT") = {5};
//+
Physical Line("FRAC") = {6};
//+
Physical Line("BC") = {1,2,3,4,5};
//+
Line Loop(1) = {1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Line{6} In Surface{1};
//+
Physical Surface("DARCY") = {1};
