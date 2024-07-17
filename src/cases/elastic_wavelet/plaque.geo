//SetFactory("OpenCASCADE");

// mesh size
h = 6e-3;

// Number of holes
N_bottom = 13;
N_side = 7;
N_top = 12;

//Radius
Radius = 0.0035;

Point(1) = {0,0,-0.0015,h};
Point(2) = {0,0.13,-0.0015,h};
Point(3) = {-0.23,0.13,-0.0015,h};
Point(4) = {-0.265,0,-0.0015,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};
//Plane Surface(1) = {1};

For r In {1:N_bottom}
    Point(5*r) = {-0.245+(r-1)*0.24/(N_bottom-1),0.005,-0.0015,h};
    Point(5*r+1) = {-0.245+(r-1)*0.24/(N_bottom-1),0.005-Radius,-0.0015,h};
    Point(5*r+2) = {-0.245+(r-1)*0.24/(N_bottom-1),0.005+Radius,-0.0015,h};
    Point(5*r+3) = {-0.245-Radius+(r-1)*0.24/(N_bottom-1),0.005,-0.0015,h};
    Point(5*r+4) = {-0.245+Radius+(r-1)*0.24/(N_bottom-1),0.005,-0.0015,h};
    Circle(4*(r-1)+5) = {5*r+1,5*r,5*r+3};
    Circle(4*(r-1)+6) = {5*r+3,5*r,5*r+2};
    Circle(4*(r-1)+7) = {5*r+2,5*r,5*r+4};
    Circle(4*(r-1)+8) = {5*r+4,5*r,5*r+1};
    Line Loop(r+1) = {4*(r-1)+5:4*(r-1)+8};
EndFor

For r In {2:N_side}
    Point(5*(r+N_bottom-1)) = {-0.005,0.005+(r-1)*0.12/(N_side-1),-0.0015,h};
    Point(5*(r+N_bottom-1)+1) = {-0.005+Radius,0.005+(r-1)*0.12/(N_side-1),-0.0015,h};
    Point(5*(r+N_bottom-1)+2) = {-0.005-Radius,0.005+(r-1)*0.12/(N_side-1),-0.0015,h};
    Point(5*(r+N_bottom-1)+3) = {-0.005,0.005+Radius+(r-1)*0.12/(N_side-1),-0.0015,h};
    Point(5*(r+N_bottom-1)+4) = {-0.005,0.005-Radius+(r-1)*0.12/(N_side-1),-0.0015,h};
    Circle(4*(r+N_bottom-1)+1) = {5*(r+N_bottom-1)+1,5*(r+N_bottom-1),5*(r+N_bottom-1)+3};
    Circle(4*(r+N_bottom-1)+2) = {5*(r+N_bottom-1)+3,5*(r+N_bottom-1),5*(r+N_bottom-1)+2};
    Circle(4*(r+N_bottom-1)+3) = {5*(r+N_bottom-1)+2,5*(r+N_bottom-1),5*(r+N_bottom-1)+4};
    Circle(4*(r+N_bottom-1)+4) = {5*(r+N_bottom-1)+4,5*(r+N_bottom-1),5*(r+N_bottom-1)+1};
    Line Loop(r+N_bottom) = {4*(r+N_bottom-1)+1:4*(r+N_bottom-1)+4};
EndFor

For r In {1:N_top-1}
    Point(5*(r+N_bottom+N_side-1)) = {-0.225+(r-1)*0.220/(N_top-1),0.125,-0.0015,h};
    Point(5*(r+N_bottom+N_side-1)+1) = {-0.225+(r-1)*0.220/(N_top-1),0.125-Radius,-0.0015,h};
    Point(5*(r+N_bottom+N_side-1)+2) = {-0.225+(r-1)*0.220/(N_top-1),0.125+Radius,-0.0015,h};
    Point(5*(r+N_bottom+N_side-1)+3) = {-0.225-Radius+(r-1)*0.220/(N_top-1),0.125,-0.0015,h};
    Point(5*(r+N_bottom+N_side-1)+4) = {-0.225+Radius+(r-1)*0.220/(N_top-1),0.125,-0.0015,h};
    Circle(4*(r+N_bottom+N_side-1)+1) = {5*(r+N_bottom+N_side-1)+1,5*(r+N_bottom+N_side-1),5*(r+N_bottom+N_side-1)+3};
    Circle(4*(r+N_bottom+N_side-1)+2) = {5*(r+N_bottom+N_side-1)+3,5*(r+N_bottom+N_side-1),5*(r+N_bottom+N_side-1)+2};
    Circle(4*(r+N_bottom+N_side-1)+3) = {5*(r+N_bottom+N_side-1)+2,5*(r+N_bottom+N_side-1),5*(r+N_bottom+N_side-1)+4};
    Circle(4*(r+N_bottom+N_side-1)+4) = {5*(r+N_bottom+N_side-1)+4,5*(r+N_bottom+N_side-1),5*(r+N_bottom+N_side-1)+1};
    Line Loop(r+N_bottom+N_side) = {4*(r+N_bottom+N_side-1)+1:4*(r+N_bottom+N_side-1)+4};
EndFor

//+
Plane Surface(2) = {1:N_bottom+N_side+N_top-1};
//+
Extrude {0, 0, 0.003} {
  Surface{2}; 
}

Point(1000) = {-0.245, 0.025, -0.0015, h};
//+
Point(1001) = {-0.125, 0.105, -0.0015, h};
//+
Point(1002) = {-0.025, 0.045, -0.0015, h};

Physical Point("Dirac1") = {1000};
Physical Point("Dirac2") = {1001};
Physical Point("Dirac3") = {1002};

Point(1003) = {-0.045, 0.085, -0.0015, 1.0};
//+
Point(1004) = {-0.185, 0.045, -0.0015, 1.0};
//+
Physical Point("PointsOfInterest1") = {1003};

Physical Point("PointsOfInterest2") = {1004};

Physical Volume("Steel") = {1};

SetElementsOrder=3;