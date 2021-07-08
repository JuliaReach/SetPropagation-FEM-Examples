l = 0.1 ;
ly = 0.4*l;

ms = 0.05*l ;
//ms = 0.5*l ;

r1  = 0.4*ly ;
x01 =  0.25*l ;
y01 =  0.5*ly ;

Point(1) = {0, 0,  0, ms};
Point(2) = {l, 0,  0, ms};
Point(3) = {l, ly, 0, ms};
Point(4) = {0.7*l, ly, 0, 0.4*ms};
Point(5) = {0.3*l, ly, 0, 0.4*ms};
Point(6) = {0, ly, 0, ms};

Point(11) = {x01-r1, y01, 0, 0.5*ms};
Point(12) = {x01   , y01, 0, 0.5*ms};
Point(13) = {x01+r1, y01, 0, 0.5*ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Circle(7) = {11, 12, 13};
Circle(8) = {13, 12, 11};

Line Loop(1) = {1, 2, 3, 4, 5, 6, -7,-8};
Line Loop(2) = {7, 8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Point ("00_01_01_00") = {3,6};
Physical Line ("00_02_03_00") = {4};

Physical Surface("01_03_00_00") = {1};
Physical Surface("02_03_00_00") = {2};
