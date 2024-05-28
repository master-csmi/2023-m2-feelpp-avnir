SetFactory("OpenCASCADE");

// mesh size
h = 0.03;

Point(1) = {0.0,0.0,0.0,h};
Point(2) = {1.0,0.0,0.0,h};
Point(3) = {0.0,-1.0,0.0,h};//+
Point(4) = {-1.0,0.0,0.0,h};//+
Point(5) = {0.0,1.0,0.0,h};//+
Circle(1) = {3, 1, 2};//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Curve Loop(1) = {4, -1, 2, 3};
//+
Curve Loop(2) = {4, -1, 2, 3};
//+
Plane Surface(1) = {2};
//+
Physical Surface("Omega")={1};
//+
Physical Curve("ABC") = {4, 1, 2, 3};
