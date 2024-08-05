SetFactory("OpenCASCADE");

// mesh size
h = 1e-1;

// square points
Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {2.0, 0.0, 0.0, h};
Point(3) = {2.0, 2.0, 0.0, h};
Point(4) = {0.0, 2.0, 0.0, h};

// square lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// square curve loop
Curve Loop(1) = {1, 2, 3, 4};

// square surface
Plane Surface(1) = {1};

// square physical groups
Physical Curve("ABC") = {1,3};
Physical Curve("Gamma") = {2,4};
Physical Surface("Omega") = {1};