// ==========================================
// Maillage du domaine Omega = (0,1) x (0,1)
// ==========================================

// Taille du Maillage
h = 0.0025;


// Sommets du carré
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

// Arêtes du carré
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};




// Boucle fermée et surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

