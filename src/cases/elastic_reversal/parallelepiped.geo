//SetFactory("OpenCASCADE");

// mesh size
h = 1e-1;

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

// square lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 6};

// square curve loop
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
// square surface
Plane Surface(1) = {1};

Extrude {0, 0, 0.01} {
  Surface{1}; 
}
Physical Curve(31) = {2, 20, 11, 16};
//+
Physical Surface("fixed") = {21};
Mesh.ElementOrder=3;//+

Physical Volume("Steel") = {1};
//+
Transfinite Surface {30} = {19, 15, 11, 10};
//+
Transfinite Surface {1} = {4, 3, 2, 1};
//+
Transfinite Surface {17} = {10, 11, 2, 1};
//+
Transfinite Surface {29} = {10, 19, 4, 1};
//+
Transfinite Surface {25} = {19, 15, 3, 4};
//+
Transfinite Surface {21} = {15, 3, 2, 11};
//+
Transfinite Curve {10, 1, 13, 4, 12, 3, 11, 2} = 200 Using Progression 1;
//+
Transfinite Curve {24, 20, 16, 15} = 1 Using Progression 1;
//+
Recombine Surface {30, 1, 17, 21, 25, 29};
