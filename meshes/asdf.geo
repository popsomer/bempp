cl__1 = 0.4;
Point(1) = {0, 0, 0.05, 0.4};
Point(2) = {0.1, 0, 0.05, 0.4};
Point(3) = {0, 0.1, 0.05, 0.4};
Point(4) = {-0.1, 0, 0.05, 0.4};
Point(5) = {0, -0.1, 0.05, 0.4};
Point(6) = {0.1, 0, 0.55, 0.4};
Point(7) = {0, 0, 0.55, 0.4};
Point(8) = {0, 0.1, 0.55, 0.4};
Point(13) = {-0.1, 0, 0.55, 0.4};
Point(18) = {0, -0.1, 0.55, 0.4};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(8) = {6, 7, 8};
Circle(9) = {8, 7, 13};
Circle(10) = {13, 7, 18};
Circle(11) = {18, 7, 6};
Line(13) = {2, 6};
Line(14) = {3, 8};
Line(18) = {4, 13};
Line(22) = {5, 18};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
Line Loop(15) = {1, 14, -8, -13};
Ruled Surface(15) = {15};
Line Loop(19) = {2, 18, -9, -14};
Ruled Surface(19) = {19};
Line Loop(23) = {3, 22, -10, -18};
Ruled Surface(23) = {23};
Line Loop(27) = {4, 13, -11, -22};
Ruled Surface(27) = {27};
Line Loop(28) = {8, 9, 10, 11};
Plane Surface(28) = {28};
Surface Loop(1) = {6, 28, 15, 19, 23, 27};
Volume(1) = {1};