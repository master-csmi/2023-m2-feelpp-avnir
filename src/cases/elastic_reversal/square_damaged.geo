//SetFactory("OpenCASCADE");

// mesh size
h = 0.03;

// square points
Point(1) = {-1.0, -1.0, 0.0, h};
Point(2) = {1.0, -1.0, 0.0, h};
Point(3) = {1.0, 1.0, 0.0, h};
Point(4) = {-1.0, 1.0, 0.0, h};
Point(5) = {-1.0, 0.0, 0.0, h};

// create a line inside the square
Point(6) = {0.3, 0.3, 0.0, h};
Point(7) = {0.3, 0.7, 0.0, h};
Point(8) = {0.4, 0.7, 0.0, h};
Point(9) = {0.4, 0.3, 0.0, h};

// create a sensor
Point(10) = {-0.9,0.,0.0,h};

// square lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 6};

Point(11) = {-0.5,0.0,0.0,h};
Point(12) = {-0.4,0.0,0.0,h};
Point(13) = {-0.6,0.0,0.0,h};

Circle(9) = {12,11,13};
Circle(10) = {13,11,12};

// square curve loop
Curve Loop(1) = {1, 2, 3, 4};
//Curve Loop(2) = {5, 6, 7, 8};
Curve Loop(2) = {9,10};
// square surface
Plane Surface(1) = {1,-2};

Physical Curve ("load") = {4};

// square physical groups
Physical Point("S1") = {10};
Physical Curve("fixed") = {2};
Physical Surface("Steel") = {1};
Mesh.ElementOrder=3;