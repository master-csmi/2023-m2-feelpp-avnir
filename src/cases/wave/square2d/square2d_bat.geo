SetFactory("OpenCASCADE");

// mesh size
h = 0.5;

// square points
Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {15.0, 0.0, 0.0, h};
Point(3) = {15.0, 10.0, 0.0, h};
Point(4) = {25.0, 10.0, 0.0, h};
Point(5) = {25.0, 0.0, 0.0, h};
Point(6) = {40.0, 0.0, 0.0, h};
Point(7) = {40.0, 20.0, 0.0, h};
Point(8) = {0.0, 20.0, 0.0, h};
Point(9) = {5.0, 5.0, 0.0, h};

// square lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

// square curve loop
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};

// square surface
Plane Surface(1) = {1};

// square physical groups
Physical Curve("Gamma") = {7};
Physical Curve("Neumann") = {7};
Physical Curve("PML") = {1,5,6,8};
Physical Surface("Omega") = {1};
Physical Point("Dirac") = {9};