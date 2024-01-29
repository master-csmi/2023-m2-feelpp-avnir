// Gmsh project created on Tue May 23 03:28:00 2023
SetFactory("OpenCASCADE");

h = 0.1; // Taille du maillage
//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {1, 0, 0, h};
//+
Point(3) = {1, 1, 0, h};
//+
Point(4) = {0, 1, 0, h};
//+
Point(5) = {0, 1, 1, h};
//+
Point(6) = {0, 0, 1, h};
//+
Point(7) = {1, 0, 1, h};
//+
Point(8) = {1, 1, 1, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {1, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 2};
//+
Line(5) = {2, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 3};
//+
Line(8) = {1, 6};
//+
Line(9) = {6, 7};
//+
Line(10) = {6, 5};
//+
Line(11) = {5, 4};
//+
Line(12) = {8, 5};
//+
Line Loop(13) = {9, 6, 12, -10};
Plane Surface(14) = {13};
//+
Line Loop(15) = {9, -5, -1, 8};
Plane Surface(16) = {15};
//+
Line Loop(17) = {2, -11, -10, -8};
Plane Surface(18) = {17};
//+
Line Loop(19) = {12, 11, 3, -7};
Plane Surface(20) = {19};
//+
Line Loop(21) = {5, 6, 7, 4};
Plane Surface(22) = {21};
//+
Line Loop(23) = {3, 4, -1, 2};
Plane Surface(24) = {23};
//+
Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25};
//+
Physical Surface("Bottom") = {24};
//+
Physical Surface("Top") = {14};
//+
Physical Surface("Front") = {20};
//+
Physical Surface("Right") = {18};
//+
Physical Surface("Back") = {16};
//+
Physical Surface("Left") = {22};
//+

// Définir le volume physique
Physical Volume("Omega") = {26};
