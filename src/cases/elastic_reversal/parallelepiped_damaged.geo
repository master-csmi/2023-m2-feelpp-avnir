//SetFactory("OpenCASCADE");
SetFactory("Built-in");

// mesh size
h = 1e-1;

// square points
Point(1) = {-1.0, -1.0, 0.0, h};
Point(2) = {1.0, -1.0, 0.0, h};
Point(3) = {1.0, 1.0, 0.0, h};
Point(4) = {-1.0, 1.0, 0.0, h};
Point(5) = {-1.0, 0.0, 0.0, h};

// square lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// square curve loop
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};


// create a cylinder
r=-0.1;
translatex = -0.6;
p5=newp; Point(p5) = { r+translatex,  0, 0, h};
p6=newp; Point(p6) = { 0+translatex,  r, 0, h};
p7=newp; Point(p7) = {-r+translatex,  0, 0, h};
p8=newp; Point(p8) = { 0+translatex, -r, 0, h};

O=newp; Point(O) = {0+translatex,0,0,h};

c5=newl; Circle(c5) = {p5,O,p6};
c6=newl; Circle(c6) = {p6,O,p7};
c7=newl; Circle(c7) = {p7,O,p8};
c8=newl; Circle(c8) = {p8,O,p5};

Line Loop(2) = {c5:c8};

Plane Surface(2) = {1, -2};

out[] = Extrude {0, 0, 0.01} {Surface{2};};
Physical Surface("fixed") = {25};
Mesh.ElementOrder=3;//+

Physical Volume("Steel") = {out[1]};
