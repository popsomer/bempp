
q = 0.4;
endpl = 0.5;
radpl = 0.1;
hwwing = 0.05;
angl = 0.3;
wings = 0.3;

Point(1) = {0,0,hwwing,q};
Point(2) = {radpl,0,hwwing,q};
Point(3) = {0,radpl,hwwing,q};
Point(4) = {-radpl,0,hwwing,q};
Point(5) = {0,-radpl,hwwing,q};
Point(6) = {radpl,0,hwwing+endpl,q};
Point(7) = {0,radpl,hwwing+endpl,q};
Point(8) = {-radpl,0,hwwing+endpl,q};
Point(9) = {0,-radpl,hwwing+endpl,q};
Point(10) = {0,0,hwwing+endpl,q};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Line(5) = {2,6};
Line(6) = {3,7};
Line(7) = {4,8};
Line(8) = {5,9};
Circle(9) = {7,10,6};
Circle(10) = {8,10,7};
Circle(11) = {9,10,8};
Circle(12) = {6,10,9};
Line(13) = {6,2};
Line(14) = {7,3};
Line(15) = {8,4};
Line(16) = {9,5};
Line Loop(17) = {1,6,9,13};
Ruled Surface(17) = {17};
Line Loop(18) = {2,7,10,14};
Ruled Surface(18) = {18};
Line Loop(19) = {3,8,11,15};
Ruled Surface(19) = {19};
Line Loop(20) = {4,5,12,16};
Ruled Surface(20) = {20};
// Above front and below back part of plane
//Translate {0, 0, -2*hwwing-endpl} { Duplicata{ Surface{17,18,19,20}; } }

// Now half-spheres for ends

Point(11) = {0,0,hwwing+radpl+endpl,q};
Circle(21) = {11,10,6};
Circle(22) = {11,10,7};
Circle(23) = {11,10,8};
Circle(24) = {11,10,9};
Circle(25) = {6,10,11};
Circle(26) = {7,10,11};
Circle(27) = {8,10,11};
Circle(28) = {9,10,11};
Line Loop(29) = {22,9,25};
Ruled Surface(29) = {29};
Line Loop(30) = {23,10,26};
Ruled Surface(30) = {30};
Line Loop(31) = {24,11,27};
Ruled Surface(31) = {31};
Line Loop(32) = {21,12,28};
Ruled Surface(32) = {32};

//Translate {0, 0, -2*hwwing-endpl} { Duplicata{ Surface{17,18,19,20}; } }

Point(12) = {radpl*Cos(angl), radpl*Sin(angl), hwwing, q};
Point(13) = {-radpl*Cos(angl), radpl*Sin(angl), hwwing, q};
Point(14) = {radpl*Cos(angl), radpl*Sin(angl), -hwwing, q};
Point(15) = {-radpl*Cos(angl), radpl*Sin(angl), -hwwing, q};
Point(16) = {0,0, -hwwing, q};
Circle(35) = {12,1,13};
//Circle(36) = {3,1,13};
Circle(36) = {15,16,14};
Line(37) = {14,12};
Line(38) = {13,15};
Line Loop(39) = {35, 38,36, 37};
Ruled Surface(39) = {39};
// Above and below parts around wings
Point(17) = {radpl*Cos(angl), -radpl*Sin(angl), hwwing, q};
Point(18) = {-radpl*Cos(angl), -radpl*Sin(angl), hwwing, q};
Point(19) = {radpl*Cos(angl), -radpl*Sin(angl), -hwwing, q};
Point(20) = {-radpl*Cos(angl), -radpl*Sin(angl), -hwwing, q};
//Point() = {0,0, -hwwing, q};
Circle(40) = {17,1,18};
//Circle(41) = {19,20,18};
Circle(41) = {20,16,19};
Line(42) = {19,17};
Line(43) = {18,20};
Line Loop(44) = {40, 43,41, 42};
Ruled Surface(44) = {44};

// below left wing

Point(21) = {radpl*Cos(angl)+wings, radpl*Sin(angl), hwwing, q};
Point(22) = {radpl*Cos(angl)+wings, -radpl*Sin(angl), hwwing, q};
Point(23) = {radpl*Cos(angl)+wings, radpl*Sin(angl), -hwwing, q};
Point(24) = {radpl*Cos(angl)+wings, -radpl*Sin(angl), -hwwing, q};
Circle(45) = {12,1,17};
Circle(46) = {14,16,19};
Line(47) = {17,22};
Line(48) = {22,21};
Line(49) = {21,12};

Line(50) = {12,21};
Line(51) = {21,23};
Line(52) = {23,14};
Line(56) = {14,12};

Line(53) = {19,24};
Line(54) = {24,23};
Line(55) = {23,14};
Line(57) = {14,19}; 

Line(58) = {19,17};
Line(60) = {22,24};
Line(61) = {24,19};

Line(63) = {23,24};
Line(64) = {24,22};

Line Loop(70) = {45,47,48,49};
Plane Surface(70) = {70};
Line Loop(71) = {50,51,52,56};
Plane Surface(71) = {71};
Line Loop(72) = {53,54,55,57};
Plane Surface(72) = {72};
Line Loop(73) = {58,47,60,61};
Plane Surface(73) = {73};
Line Loop(74) = {51,63,64,48};
Plane Surface(74) = {74};


Symmetry{ 0.0, 0.0, 1.0, 0.0 }{Duplicata{Surface{17,18,19,20,29,30,31,32};}}


Symmetry{ 1.0, 0.0, 0.0, 0.0 }{Duplicata{Surface{70,71,72,73,74};}}


/*
give error
Point(11) = {0,0,hwwing+radpl,q};
Circle(21) = {11,1,2};
Circle(22) = {11,1,3};
Circle(23) = {11,1,4};
Circle(24) = {11,1,5};
Circle(25) = {2,1,11};
Circle(26) = {3,1,11};
Circle(27) = {4,1,11};
Circle(28) = {5,1,11};
Line Loop(29) = {21,1,26};
Ruled Surface(29) = {29};
Line Loop(30) = {22,2,27};
Ruled Surface(30) = {30};
Line Loop(31) = {23,3,28};
Ruled Surface(31) = {31};
Line Loop(32) = {24,4,25};
Ruled Surface(32) = {32};
*/
//Translate {0,0,endpl} {Surface{29,30,31,32}; }

//Translate {0, 0, -2*hwwing-endpl} { Duplicata{ Surface{17,18,19,20}; } }


//Line Loop(5) = {1,2,3,4};
//Plane Surface(6) = {5};
//Plane Surface(6) = {1,2,3,4};
//Extrude {0,0,endpl} {
//  Surface{6};
//  Surface{5};
//}
//Delete{Surface{6}; }

/*
// Above front part plane, below back
Point(6) = {0,0,-hwwing,q};
Point(7) = {radpl,0,-hwwing,q};
Point(8) = {0,radpl,-hwwing,q};
Point(9) = {-radpl,0,-hwwing,q};
Point(10) = {0,-radpl,-hwwing,q};
Circle(7) = {7,6,8};
Circle(8) = {8,6,9};
Circle(9) = {9,6,10};
Circle(10) = {10,6,7};
Line Loop(11) = {7,8,9,10};
Plane Surface(12) = {11};
Extrude {0,0,-endpl} {
  Surface{12};
}
*/
// Now parts around wings
//Point(11) = {radpl*cos(angl), radpl*sin(angl),-hwwing,q};
//Point(12) = {radpl*cos(angl), radpl*sin(angl),-hwwing,q};
//Circle(13) = {7,6,8};
//Circle(14) = {8,6,9};
//Line Loop(11) = {7,8,9,10};
//Plane Surface(12) = {11};
//Extrude {0,0,-endpl} {
//  Surface{12};
//}
