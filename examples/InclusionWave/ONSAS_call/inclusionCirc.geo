l = 0.1 ;
ly = 0.5*l;

ms = 0.05*l ;
//ms = 0.5*l ;

r1  = 0.2*ly ;
x01 =  0.5*l-r1 ;
y01 =  0.5*ly ;

Point(1) = {0, 0,  0, ms};
Point(2) = {l, 0,  0, ms};
Point(3) = {l, ly, 0, ms};
Point(4) = {0.5*l, ly, 0, 0.5*ms};
Point(5) = {0, ly, 0, ms};

Point(6) = {x01-r1, y01, 0, 0.5*ms};
Point(7) = {x01   , y01, 0, 0.5*ms};
Point(8) = {x01+r1, y01, 0, 0.5*ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

Circle(6) = {6, 7, 8};
Circle(7) = {8, 7, 6};

Line Loop(1) = {1, 2, 3, 4,5,-6,-7};
Line Loop(2) = {6, 7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Point ("00_01_01_00") = {1};
Physical Point ("00_01_02_00") = {2};
Physical Point ("00_01_03_00") = {4};

Physical Surface("01_03_00_00") = {1};
Physical Surface("02_03_00_00") = {2};
