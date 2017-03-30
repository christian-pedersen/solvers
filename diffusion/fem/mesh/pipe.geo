lco = 0.05;
lci = 0.1;
x0 = 0;
y0 = 0;


//center
Point(1) = {x0, y0, 0, lc};

//outer edge
ro = 1.0;
Point(2) = {x0+ro, y0, 0, lco};
Point(3) = {x0, yc0+ro, 0, lco};
Point(4) = {x0-ro, y0, 0, lco};
Point(5) = {x0, y0-ro, 0, lco};

//inner edge
ri = 0.4;
Point(6) = {x0+ri, y0, 0, lci};
Point(7) = {x0, yc0+ri, 0, lci};
Point(8) = {x0-ri, y0, 0, lci};
Point(9) = {x0, y0-ri, 0, lci};


// Lines for outer edge
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Line Loop(5) = {1, 2, 3, 4};

// Lines for inner edge
Circle(6) = {6, 1, 7};
Circle(7) = {7, 1, 8};
Circle(8) = {8, 1, 9};
Circle(9) = {9, 1, 6};

Line Loop(10) = {6, 7, 8, 9};

Plane Surface(11) = {5, 10};
Plane Surface(12) = {10};
