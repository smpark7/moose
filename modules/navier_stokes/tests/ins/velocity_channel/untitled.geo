//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 75, 0, 1.0};
//+
Point(3) = {0, -75, 0, 1.0};
//+
Point(4) = {200, -75, 0, 1.0};
//+
Point(5) = {200, 0, 0, 1.0};
//+
Point(6) = {200, 75, 0, 1.0};
//+ Top
Line(1) = {2, 6};
//+ top half of right boundary
Line(2) = {6, 5};
//+ bottom half of right boundary
Line(3) = {5, 4};
//+ bottom
Line(4) = {4, 3};
//+ bottom half of left boundary
Line(5) = {3, 1};
//+ top half of left boundary
Line(6) = {1, 2};
//+
Line Loop(7) = {1, 2, 3, 4, 5, 6};
//+
Plane Surface(8) = {7};

Transfinite Line{1, 4} = 31;
Transfinite Line{2, 3, 5, 6} = 16;
Transfinite Surface{8};
Recombine Surface{8};
