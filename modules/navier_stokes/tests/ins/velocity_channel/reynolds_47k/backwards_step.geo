element_size = 0.5;

// Top right
Point(3) = {24, 0, 0, 1.0};
// Top step
Point(4) = {4, 0, 0, 1.0};
// Step left
Point(5) = {0, -2, 0, 1.0};
// Top left
Point(6) = {0, 0, 0, 1.0};
// Step step
Point(7) = {4, -2, 0, 1.0};
// Bottom step
Point(8) = {4, -3, 0, 1.0};
// Bottom right
Point(9) = {24, -3, 0, 1.0};
// Step right
Point(10) = {24, -2, 0, 1.0};

// LHS of top boundary
Line(1) = {6, 4};
// RHS of top boundary
Line(2) = {4, 3};
// Upper half of right boundary
Line(3) = {3, 10};
// Bottom half of right boundary
Line(4) = {10, 9};
// Bottom of expanded channel
Line(5) = {9, 8};
// Step
Line(6) = {8, 7};
// Bottom of inlet channel
Line(7) = {7, 5};
// Left boundary / inlet
Line(8) = {5, 6};
// End of inlet channel
Line(9) = {4, 7};
// Divider of original channel and expanded channel
Line(10) = {7, 10};

Line Loop(11) = {1, 9, 7, 8};
Plane Surface(12) = {11};
Line Loop(13) = {2, 3, -10, -9};
Plane Surface(14) = {13};
Line Loop(15) = {5, 6, 10, 4};
Plane Surface(16) = {15};

Transfinite Line{1, 7} = 4. / element_size + 1;
Transfinite Line{2, 5, 10} = 20. / element_size + 1;
Transfinite Line{3, 8, 9} = 2. / element_size + 1;
Transfinite Line{4, 6} = 1. / element_size + 1;
Transfinite Surface{12,14,16};
Recombine Surface{12,14,16};

Physical Point("corner") = {3};

Physical Line("inlet") = {8};
Physical Line("norm_y_walls") = {1, 2, 5, 7};
Physical Line("norm_x_walls") = {6};
Physical Line("outlet") = {3, 4};

Physical Surface("domain") = {12, 14, 16};
