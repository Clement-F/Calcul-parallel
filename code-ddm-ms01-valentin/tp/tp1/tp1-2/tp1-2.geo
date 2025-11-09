// ===================================================
// Maillage du domaine Omega = (0,1) x (0,1)x(0,1)
// ===================================================


// Taille du Maillage
h = 0.1;

// Sommets du cube
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {0, 0, 1, h};
Point(6) = {1, 0, 1, h};
Point(7) = {1, 1, 1, h};
Point(8) = {0, 1, 1, h};


// Arêtes du cube
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};


// Faces du cube
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};

Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {15};

Line Loop(17) = {1, 10, -5, -9};
Plane Surface(18) = {17};

Line Loop(19) = {2, 11, -6, -10};
Plane Surface(20) = {19};

Line Loop(21) = {3, 12, -7, -11};
Plane Surface(22) = {21};

Line Loop(23) = {4, 9, -8, -12};
Plane Surface(24) = {23};


// Volume
Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25};

Mesh 3;